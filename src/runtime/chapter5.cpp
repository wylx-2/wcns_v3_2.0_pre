#include <wcns/runtime/chapter5.hpp>
#include <wcns/mesh/linear_operators.hpp>
#include <wcns/solver/viscous_boundary.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <functional>
#include <limits>
#include <sys/stat.h>

namespace wcns {
namespace {
constexpr int NG=3,NW=8,NM=24,WM=13;
const PartitionLeaf& leaf(const StatisticContext& c,BlockId id) {
    for(const auto& l:c.partition.leaves()) if(l.block==id) return l;
    throw std::runtime_error("chapter5: missing leaf");
}
void root_io(const MpiRuntime& mpi,const std::function<void()>& action) {
    std::string error;
    if(mpi.rank()==0) try { action(); } catch(const std::exception& e) {error=e.what();}
    error=mpi.broadcast_string(std::move(error));
    if(!error.empty()) throw std::runtime_error(error);
}
std::ofstream file(const std::string& path,bool append=false) {
    std::ofstream f; f.exceptions(std::ios::failbit|std::ios::badbit);
    f.open(path,append?std::ios::app:std::ios::trunc); f<<std::setprecision(17);return f;
}
bool exists(const std::string& p) {struct stat s{};return stat(p.c_str(),&s)==0;}
Real variance(Real raw,Real mean) {return std::max(Real(0),raw-mean*mean);}
}
const Chapter5Runtime::Zone& Chapter5Runtime::zone(BlockId id) const {
    for(const auto& z:zones_) if(z.id==id) return z;
    throw std::runtime_error("chapter5: missing zone");
}
Chapter5Runtime::Chapter5Runtime(const CaseConfig& cfg,const StatisticContext& c):config_(cfg),context_(c) {
    std::size_t n=0,nw=0;
    for(const auto& z:c.partition.zones()) {
        const auto e=z.cell_extent;
        if(z.cell_dimension!=3 || e.nj<3) throw std::invalid_argument("chapter5 requires 3D extruded zones");
        zones_.push_back({z.source_zone,e.ni,e.nj,e.nk,n,nw});n+=e.ni*e.nj;nw+=e.ni;
    }
    geometry_.assign(n*NG,0);walls_.assign(nw*NW,0);
    instant_.assign(n*NM,0);integral_=instant_;wall_.assign(nw*WM,0);wall_integral_=wall_;
    bool valid=true;
    const Real span=cfg.chapter5.type=="sd7003"?.2:6.;
    for(const auto& b:c.local_blocks.blocks()) {
        const auto& l=leaf(c,b.id());const auto& z=zone(l.source_zone);const auto e=b.cell_extent();
        const auto& p=b.coordinates;const auto& cc=c.metrics.at(b.id()).cell_coordinates();
        const BoundaryPatch* wall_patch=nullptr;
        if(l.cells.begin.j==0) {
            valid=valid && e.nj>=(c.profile.kind()==AlgorithmProfileKind::PhengleiWcns?4:6);
            bool wall=false;
            for(const auto& patch:b.boundaries) if(patch.face.axis==Axis::J && patch.face.side==Side::Lower) {
                wall_patch=&patch;
                wall=patch.type==(cfg.chapter5.type=="compression_ramp"?BoundaryType::NoSlipIsothermalWall:BoundaryType::NoSlipAdiabaticWall);
                if(c.boundary_data) {
                    const auto& data=c.boundary_data->at(b.id()).at(patch.name);
                    for(Real u:data.wall_velocity) valid=valid && u==0;
                    if(cfg.chapter5.type=="compression_ramp") valid=valid && data.wall_temperature
                        && std::abs(*data.wall_temperature-323./170.)<1e-12;
                }
            }
            valid=valid && wall;
        }
        for(int k=0;k<e.nk;++k) for(int j=0;j<e.nj;++j) for(int i=0;i<e.ni;++i) {
            const auto q=z.offset+(j+l.cells.begin.j)*z.nx+i+l.cells.begin.i;
            const Real dz=p.z(i,j,k+1)-p.z(i,j,k);
            valid=valid && dz>0 && std::abs(dz-span/z.nz)<1e-10
                && std::abs(p.x(i,j,k+1)-p.x(i,j,k))<1e-10 && std::abs(p.y(i,j,k+1)-p.y(i,j,k))<1e-10;
            geometry_[q*NG]+=cc.x(i,j,k)*dz;geometry_[q*NG+1]+=cc.y(i,j,k)*dz;geometry_[q*NG+2]+=dz;
            if(j+l.cells.begin.j==0) {
                auto w=(z.wall_offset+i+l.cells.begin.i)*NW;
                const auto& area=c.metrics.at(b.id()).j_faces();
                const Real ax=area.x(i,j,k),ay=area.y(i,j,k),length=std::hypot(ax,ay),ds=length/dz;
                valid=valid && ds>0;
                const Real nx=ax/length,ny=ay/length,sign=ny>=0?1.:-1.;
                if(!wall_patch) {valid=false;continue;}
                const auto wp=viscous_wall_coordinates(b,c.metrics.at(b.id()),*wall_patch,{i,j,k},c.profile);
                const Real values[NW]={wp[0],wp[1],sign*ny,-sign*nx,nx,ny,ds,span};
                for(int m=0;m<NW;++m) walls_[w+m]+=dz*values[m];
            }
        }
    }
    if(!c.mpi.all_true(valid)) throw std::invalid_argument("chapter5 expects uniform-z extrusion with stationary jmin walls");
    c.mpi.sum_reals(geometry_);c.mpi.sum_reals(walls_);
    for(std::size_t q=0;q<n;++q) {
        if(std::abs(geometry_[q*NG+2]-span)>1e-9) throw std::invalid_argument("chapter5 span mismatch");
        geometry_[q*NG]/=span;geometry_[q*NG+1]/=span;
    }
    for(auto& v:walls_) v/=span;
    for(std::size_t p=0;p<cfg.chapter5.probes.size();p+=3) {
        const Real x=cfg.chapter5.probes[p],y=cfg.chapter5.probes[p+1],zz=cfg.chapter5.probes[p+2];
        Real best=1e300;std::size_t cell=0;int nz=0;
        for(const auto& z:zones_) for(int j=0;j<z.ny;++j) for(int i=0;i<z.nx;++i) {
            const auto q=z.offset+j*z.nx+i;const Real d=std::hypot(x-geometry_[q*NG],y-geometry_[q*NG+1]);
            if(d<best) {best=d;cell=q;nz=z.nz;}
        }
        if(zz<0 || zz>span) throw std::invalid_argument("probe z outside span");
        const int k=std::min(nz-1,static_cast<int>(zz/span*nz));
        probes_.push_back({cell,k,geometry_[cell*NG],geometry_[cell*NG+1],span*(k+.5)/nz});
    }
    probe_values_.assign(probes_.size()*6,0);
}
void Chapter5Runtime::sample() {
    const auto& c=context_;const auto& ref=c.quantities.reference;const auto& gas=c.quantities.gas;
    std::fill(instant_.begin(),instant_.end(),0);std::fill(wall_.begin(),wall_.end(),0);std::fill(probe_values_.begin(),probe_values_.end(),0);
    bool valid=true;
    auto state=[&](const StructuredBlock& b,int i,int j,int k) {
        return temperature_primitive_from_conservative(load_conservative(b.flow.conservative,{i,j,k}),gas,ref,c.quantities.floors,3);
    };
    for(const auto& b:c.local_blocks.blocks()) {
        const auto& l=leaf(c,b.id());const auto& z=zone(l.source_zone);const auto e=b.cell_extent();
        const auto& cc=c.metrics.at(b.id()).cell_coordinates();
        for(int k=0;k<e.nk;++k) for(int j=0;j<e.nj;++j) for(int i=0;i<e.ni;++i) {
            const auto q=z.offset+(j+l.cells.begin.j)*z.nx+i+l.cells.begin.i;
            const Real dz=b.coordinates.z(i,j,k+1)-b.coordinates.z(i,j,k);
            const auto p=state(b,i,j,k);const Real rho=p[0],u=p[1],v=p[2],w=p[3],T=p[4],P=rho*T/(gas.gamma()*ref.mach()*ref.mach());
            const Real a[NM]={rho,u,v,w,P,T,u*u,v*v,w*w,u*v,u*w,v*w,rho*u,rho*v,rho*w,
                rho*u*u,rho*v*v,rho*w*w,rho*u*v,rho*u*w,rho*v*w,rho*rho,P*P,T*T};
            for(int m=0;m<NM;++m) {valid=valid && std::isfinite(a[m]);instant_[q*NM+m]+=dz*a[m];}
            for(std::size_t h=0;h<probes_.size();++h) if(probes_[h].cell==q && probes_[h].k==k+l.cells.begin.k)
                for(int m=0;m<6;++m) probe_values_[h*6+m]+=a[m];
            if(j+l.cells.begin.j!=0) continue;
            const auto iw=z.wall_offset+i+l.cells.begin.i;const auto* g=walls_.data()+iw*NW;
            const Real d1=(cc.x(i,j,k)-g[0])*g[4]+(cc.y(i,j,k)-g[1])*g[5];
            if(!(d1>0)) {valid=false;continue;}
            // Evaluate the accepted conservative state, not a stale RK-stage cache.
            // Match the solver's wall interpolation and Dirichlet closure.
            Real pw=0,Tface=0;
            const auto& row=cached_line_operators(c.profile,e.nj).interpolation_rows()[0];
            for(const auto [jj,coefficient]:row) {
                const auto pwall=state(b,i,jj,k);
                pw+=coefficient*pwall[0]*pwall[4]/(gas.gamma()*ref.mach()*ref.mach());
                Tface+=coefficient*pwall[4];
            }
            const bool ramp=config_.chapter5.type=="compression_ramp";
            const Real Tw=ramp?323./170.:Tface;
            if(!(Tw>0 && pw>0)) {valid=false;continue;}
            const int count=c.profile.kind()==AlgorithmProfileKind::PhengleiWcns?4:6;
            std::vector<Real> tangent(count),spanwise(count),temperature(count);
            for(int jj=0;jj<count;++jj) {
                const auto s=state(b,i,jj,k);
                tangent[jj]=s[1]*g[2]+s[2]*g[3];spanwise[jj]=s[3];temperature[jj]=s[4];
            }
            const Real scale=.5/d1;
            const Real mu=c.quantities.transport.viscosity(Tw)/ref.reynolds();
            const Real rho_w=pw*gas.gamma()*ref.mach()*ref.mach()/Tw;
            const Real tau=mu*scale*wall_dirichlet_computational_derivative(0,tangent,c.profile);
            const Real tau_z=mu*scale*wall_dirichlet_computational_derivative(0,spanwise,c.profile);
            const Real heat=ramp?c.quantities.transport.thermal_coefficient(Tw,gas,ref)/ref.reynolds()
                *scale*wall_dirichlet_computational_derivative(Tw,temperature,c.profile):0;
            const Real ut=std::sqrt(std::hypot(tau,tau_z)/rho_w),yp=ut*d1*rho_w/mu;
            const Real vwall[WM]={pw,pw*pw,tau,tau_z,heat,heat*heat,Tw,rho_w,ut,yp,tau<0?1.:0.,mu,d1};
            for(int m=0;m<WM;++m) {valid=valid && std::isfinite(vwall[m]);wall_[iw*WM+m]+=dz*vwall[m];}
        }
    }
    if(!c.mpi.all_true(valid)) throw std::runtime_error("chapter5: invalid sample or wall geometry");
    c.mpi.sum_reals(instant_);c.mpi.sum_reals(wall_);c.mpi.sum_reals(probe_values_);
    for(std::size_t q=0;q<instant_.size()/NM;++q) for(int m=0;m<NM;++m) instant_[q*NM+m]/=geometry_[q*NG+2];
    loads_={0,0,0};const Real pinf=1/(gas.gamma()*ref.mach()*ref.mach());
    for(std::size_t q=0;q<wall_.size()/WM;++q) {
        const auto* g=walls_.data()+q*NW;auto* a=wall_.data()+q*WM;
        for(int m=0;m<WM;++m) a[m]/=g[7];
        const Real fx=(-(a[0]-pinf)*g[4]+a[2]*g[2])*g[6],fy=(-(a[0]-pinf)*g[5]+a[2]*g[3])*g[6];
        loads_[0]+=2*fx;loads_[1]+=2*fy;loads_[2]+=2*((g[0]-.25)*fy-g[1]*fx);
    }
    if(config_.chapter5.type=="sd7003") {
        const Real a=4*std::acos(-1.)/180,fx=loads_[0],fy=loads_[1];
        loads_[0]=fx*std::cos(a)+fy*std::sin(a);loads_[1]=-fx*std::sin(a)+fy*std::cos(a);
    }
}
void Chapter5Runtime::on_initial(const SimulationState& s) {
    sample();if(!restored_) {last_step_=s.step;last_time_=s.time;}
    if(last_step_!=s.step || std::abs(last_time_-s.time)>1e-12) throw std::runtime_error("chapter5 checkpoint time mismatch");
}
void Chapter5Runtime::on_step(const SimulationState& s,bool) {
    if(!s.residuals.finite || s.step<=last_step_ || s.time<=last_time_) return;
    sample();const auto& b=config_.chapter5;
    const Real dt=std::max(Real(0),std::min(s.time,b.statistics_end)-std::max(last_time_,b.statistics_start));
    if(dt>0) {
        for(std::size_t i=0;i<integral_.size();++i) integral_[i]+=dt*instant_[i];
        for(std::size_t i=0;i<wall_.size();++i) wall_integral_[i]+=dt*wall_[i];
        for(int i=0;i<3;++i) {load_integral_[i]+=dt*loads_[i];load_integral_[i+3]+=dt*loads_[i]*loads_[i];}
        weight_+=dt;++samples_;
    }
    last_step_=s.step;last_time_=s.time;
    if(s.step%b.history_every_steps==0) write_history(s);
    if(s.step%b.write_every_steps==0) write_means();
}
void Chapter5Runtime::on_final(const SimulationState& s) {
    if(s.stop_reason==StopReason::NumericalFailure) return;
    if(s.step!=last_history_) write_history(s);write_means();
}
std::vector<std::string> Chapter5Runtime::output_paths() const {
    std::vector<std::string> p;
    for(const char* s:{"history","probes","mean_xy","wall","loads"}) p.push_back(config_.output.directory+"/"+config_.case_name+"_benchmark_"+s+".csv");
    return p;
}
void Chapter5Runtime::write_history(const SimulationState& s) {
    root_io(context_.mpi,[&] {
        auto paths=output_paths();const bool append=exists(paths[0]);auto f=file(paths[0],append);
        if(!append) f<<"step,time,Cd_or_Cx,Cl_or_Cy,Cmz,weight,samples,xs,xr,wall_pressure_gradient_x,q_into_wall_max\n";
        const bool ramp=config_.chapter5.type=="compression_ramp";
        std::vector<std::size_t> selected;
        for(std::size_t q=0;q<wall_.size()/WM;++q) {
            const auto* g=walls_.data()+q*NW;
            if(ramp?(g[0]>75 && g[0]<140):(g[5]>0 && g[0]>.05 && g[0]<.98)) selected.push_back(q);
        }
        std::sort(selected.begin(),selected.end(),[&](auto a,auto b){return walls_[a*NW]<walls_[b*NW];});
        Real xs=std::numeric_limits<Real>::quiet_NaN(),xr=xs,shock=xs,qmax=-1e300;
        Real start=xs,longest=0,gradient=0;
        for(std::size_t i=0;i<selected.size();++i) {
            const auto q=selected[i];qmax=std::max(qmax,wall_[q*WM+4]);
            if(!i) continue;
            const auto p=selected[i-1];const Real x0=walls_[p*NW],x1=walls_[q*NW];
            const Real a=wall_[p*WM+2],b=wall_[q*WM+2];
            if(a>0 && b<0) start=x0-a*(x1-x0)/(b-a);
            if(a<0 && b>0 && std::isfinite(start)) {
                const Real end=x0-a*(x1-x0)/(b-a);
                if(end-start>longest) {xs=start;xr=end;longest=end-start;}
            }
            const Real dp=(wall_[q*WM]-wall_[p*WM])/(x1-x0);
            if(ramp && x1<125 && dp>gradient) {gradient=dp;shock=.5*(x0+x1);}
        }
        f<<s.step<<','<<s.time;for(Real v:loads_) f<<','<<v;
        f<<','<<weight_<<','<<samples_<<','<<xs<<','<<xr<<','<<shock<<','<<qmax<<'\n';
        const bool ap=exists(paths[1]);auto p=file(paths[1],ap);
        if(!ap) p<<"step,time,probe,x,y,z,rho,u,v,w,p,T\n";
        for(std::size_t q=0;q<probes_.size();++q) {
            const auto& a=probes_[q];p<<s.step<<','<<s.time<<','<<q<<','<<a.x<<','<<a.y<<','<<a.z;
            for(int m=0;m<6;++m) p<<','<<probe_values_[q*6+m];p<<'\n';
        }
    });last_history_=s.step;
}
void Chapter5Runtime::write_means() const {
    root_io(context_.mpi,[&] {
        auto paths=output_paths();auto f=file(paths[2]),w=file(paths[3]),l=file(paths[4]);
        f<<"zone,i,j,x,y,wall_x,wall_y,tx,ty,nx,ny,rho,u,v,w,p,T,uu,vv,ww,uv,uw,vw,favre_u,favre_v,favre_w,favre_uu,favre_vv,favre_ww,favre_uv,favre_uw,favre_vw,rho_rms,p_rms,T_rms,weight,samples\n";
        w<<"zone,i,x,y,tx,ty,nx,ny,ds,p,Cp,p_rms,tau_t,Cf_t,Cfx,tau_z,q_into_wall,q_rms,Twall,rho_wall,u_tau_mean,yplus_mean,backflow_fraction,mu_wall_over_Re,first_cell_distance,weight\n";
        l<<"coefficient,mean,rms,weight,samples\n";
        if(weight_<=0) return;
        const int c1[6]={1,2,3,1,1,2},c2[6]={1,2,3,2,3,3};
        for(const auto& z:zones_) {
            for(int j=0;j<z.ny;++j) for(int i=0;i<z.nx;++i) {
                const auto q=z.offset+j*z.nx+i,iw=z.wall_offset+i;const auto* g=walls_.data()+iw*NW;
                Real a[NM];for(int m=0;m<NM;++m) a[m]=integral_[q*NM+m]/weight_;
                f<<z.id<<','<<i<<','<<j<<','<<geometry_[q*NG]<<','<<geometry_[q*NG+1];
                for(int m=0;m<6;++m) f<<','<<g[m];for(int m=0;m<6;++m) f<<','<<a[m];
                for(int m=0;m<6;++m) f<<','<<(a[6+m]-a[c1[m]]*a[c2[m]]);
                for(int m=0;m<3;++m) f<<','<<a[12+m]/a[0];
                for(int m=0;m<6;++m) f<<','<<(a[15+m]/a[0]-a[11+c1[m]]*a[11+c2[m]]/(a[0]*a[0]));
                for(int m=0;m<3;++m) f<<','<<std::sqrt(variance(a[21+m],a[m==0?0:m==1?4:5]));
                f<<','<<weight_<<','<<samples_<<'\n';
            }
            const Real pinf=1/(context_.quantities.gas.gamma()*std::pow(context_.quantities.reference.mach(),2));
            for(int i=0;i<z.nx;++i) {
                const auto q=z.wall_offset+i;const auto* g=walls_.data()+q*NW;
                Real a[WM];for(int m=0;m<WM;++m) a[m]=wall_integral_[q*WM+m]/weight_;
                w<<z.id<<','<<i;for(int m=0;m<7;++m) w<<','<<g[m];
                w<<','<<a[0]<<','<<2*(a[0]-pinf)<<','<<std::sqrt(variance(a[1],a[0]))
                 <<','<<a[2]<<','<<2*a[2]<<','<<2*a[2]*g[2]<<','<<a[3]<<','<<a[4]<<','<<std::sqrt(variance(a[5],a[4]));
                for(int m=6;m<WM;++m) w<<','<<a[m];w<<','<<weight_<<'\n';
            }
        }
        const char* names[3]={"Cd_or_Cx","Cl_or_Cy","Cmz"};
        for(int m=0;m<3;++m) l<<names[m]<<','<<load_integral_[m]/weight_<<','<<std::sqrt(variance(load_integral_[m+3]/weight_,load_integral_[m]/weight_))<<','<<weight_<<','<<samples_<<'\n';
    });
}
std::string Chapter5Runtime::serialize() const {
    std::ostringstream s;s<<std::setprecision(17)<<config_.chapter5.signature()<<'\n'
        <<integral_.size()<<' '<<wall_integral_.size()<<' '<<weight_<<' '<<samples_<<' '<<last_step_<<' '<<last_time_<<'\n';
    for(const auto* a:{&integral_,&wall_integral_}) {for(Real v:*a) s<<v<<' ';s<<'\n';}
    for(Real v:load_integral_) s<<v<<' ';s<<'\n';return s.str();
}
void Chapter5Runtime::restore(const std::string& text) {
    std::istringstream s(text);std::string signature;std::getline(s,signature);
    std::size_t n=0,nw=0,samples=0,step=0;Real weight=0,time=0;s>>n>>nw>>weight>>samples>>step>>time;
    if(!s || signature!=config_.chapter5.signature() || n!=integral_.size() || nw!=wall_integral_.size()
       || !std::isfinite(weight) || weight<0 || !std::isfinite(time) || time<0) throw std::runtime_error("invalid chapter5 checkpoint header");
    auto a=integral_,b=wall_integral_;auto loads=load_integral_;
    for(auto* v:{&a,&b}) for(Real& x:*v) {s>>x;if(!s || !std::isfinite(x)) throw std::runtime_error("invalid chapter5 moments");}
    for(Real& x:loads) {s>>x;if(!s || !std::isfinite(x)) throw std::runtime_error("invalid chapter5 loads");}
    std::string extra;if(s>>extra) throw std::runtime_error("trailing chapter5 state");
    integral_=std::move(a);wall_integral_=std::move(b);load_integral_=loads;
    weight_=weight;samples_=samples;last_step_=step;last_time_=time;restored_=true;
}
} // namespace wcns
