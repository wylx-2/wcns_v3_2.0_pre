#include <wcns/runtime/periodic_hill.hpp>
#include <wcns/physics/periodic_hill.hpp>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace wcns {
namespace {
constexpr int ng = 7, nm = 13;
constexpr Real area = 2.035 * 4.5;
const PartitionLeaf& leaf_for(const StatisticContext& c, BlockId id) {
    for (const auto& l : c.partition.leaves()) if (l.block == id) return l;
    throw std::runtime_error("hill: missing partition leaf");
}
// Propagate root I/O errors before the next collective solver operation.
void root_io(const MpiRuntime& mpi, const std::function<void()>& action) {
    std::string error;
    if (mpi.rank() == 0) try { action(); } catch (const std::exception& e) { error=e.what(); }
    error=mpi.broadcast_string(std::move(error));
    if (!error.empty()) throw std::runtime_error(error);
}
std::ofstream output_file(const std::string& path, bool append=false) {
    std::ofstream f;
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f.open(path, append ? std::ios::app : std::ios::trunc);
    f << std::setprecision(17);
    return f;
}
}

PeriodicHillRuntime::PeriodicHillRuntime(const CaseConfig& config, const StatisticContext& c)
    : config_(config), context_(c), force_(config.source_terms.pressure_gradient[0]) {
    if (c.partition.zones().size()!=1 || c.partition.zones()[0].cell_dimension!=3)
        throw std::invalid_argument("hill requires one extruded 3D structured source zone");
    const auto ext=c.partition.zones()[0].cell_extent;
    nx_=ext.ni; ny_=ext.nj;
    if (nx_<8 || ny_<4 || ext.nk<4) throw std::invalid_argument("hill grid is too small");
    geometry_.assign(static_cast<std::size_t>(nx_)*ny_*ng,0);
    instant_.assign(static_cast<std::size_t>(nx_)*ny_*nm,0);
    integral_=instant_;
    wall_.assign(static_cast<std::size_t>(nx_)*4,0); wall_integral_=wall_;
    bool valid=true;
    Real lower_patches=0,upper_patches=0;
    for (const auto& b:c.local_blocks.blocks()) {
        const auto& l=leaf_for(c,b.id());
        const auto& nodes=b.coordinates;
        const auto& cc=c.metrics.at(b.id()).cell_coordinates();
        const auto e=b.cell_extent();
        if(c.boundary_data) {
            for(Axis axis:{Axis::I,Axis::K}) for(Side side:{Side::Lower,Side::Upper})
                valid=valid && connection_side_is_fully_covered(b,axis,side);
            for(const auto& patch:b.boundaries) {
                valid=valid && patch.face.axis==Axis::J
                    && (patch.type==BoundaryType::NoSlipIsothermalWall || patch.type==BoundaryType::NoSlipAdiabaticWall);
                if(patch.face.side==Side::Lower) ++lower_patches; else ++upper_patches;
                const auto& data=c.boundary_data->at(b.id()).at(patch.name);
                for(Real v:data.wall_velocity) valid=valid && v==0;
            }
        }
        for(int k=0;k<e.nk;++k) for(int j=0;j<e.nj;++j) for(int i=0;i<e.ni;++i) {
            const int I=i+l.cells.begin.i,J=j+l.cells.begin.j;
            const auto q=static_cast<std::size_t>(J*nx_+I)*ng;
            const Real dz=nodes.z(i,j,k+1)-nodes.z(i,j,k);
            valid=valid && dz>0 && std::abs(nodes.x(i,j,k)-nodes.x(i,j,k+1))<1e-10
                && std::abs(nodes.y(i,j,k)-nodes.y(i,j,k+1))<1e-10
                && std::abs(nodes.x(i,j,k)-nodes.x(i,j+1,k))<1e-10
                && std::abs(nodes.z(i,j,k)-nodes.z(i+1,j,k))<1e-10
                && std::abs(nodes.z(i,j,k)-nodes.z(i,j+1,k))<1e-10;
            geometry_[q]+=cc.x(i,j,k)*dz; geometry_[q+1]+=cc.y(i,j,k)*dz;
            geometry_[q+2]+=dz;
            if(J==0) {
                geometry_[q+3]+=.5*(nodes.y(i,0,k)+nodes.y(i+1,0,k))*dz;
                geometry_[q+5]+=(nodes.y(i+1,0,k)-nodes.y(i,0,k))
                    /(nodes.x(i+1,0,k)-nodes.x(i,0,k))*dz;
                valid=valid && std::abs(nodes.y(i,0,k)-periodic_hill_geometry(nodes.x(i,0,k))[0])<1e-9;
            }
            if(J==ny_-1) geometry_[q+4]+=.5*(nodes.y(i,j+1,k)+nodes.y(i+1,j+1,k))*dz;
            if(I==0) {
                valid=valid && std::abs(nodes.x(i,j,k))<1e-10;
                geometry_[q+6]+=(nodes.y(i,j+1,k)-nodes.y(i,j,k))*dz;
            }
            if(I==nx_-1) valid=valid && std::abs(nodes.x(i+1,j,k)-9)<1e-10;
        }
    }
    if(c.boundary_data) {
        const Real lower=c.mpi.sum(lower_patches),upper=c.mpi.sum(upper_patches);
        valid=valid && lower>0 && upper>0;
    }
    if(!c.mpi.all_true(valid)) throw std::invalid_argument("hill requires the canonical extruded 9 x 3.035 x 4.5 mesh");
    c.mpi.sum_reals(geometry_);
    for(int j=0;j<ny_;++j) for(int i=0;i<nx_;++i) {
        const auto q=static_cast<std::size_t>(j*nx_+i)*ng;
        const Real w=geometry_[q+2];
        if(std::abs(w-4.5)>1e-9) throw std::invalid_argument("hill span must equal 4.5");
        for(int a:{0,1,3,4,5}) geometry_[q+a]/=w;
    }
    Real crest_area=0;
    for(int j=0;j<ny_;++j) {
        crest_area+=geometry_[static_cast<std::size_t>(j*nx_)*ng+6];
        for(int i=0;i<nx_;++i) {
            const auto q=static_cast<std::size_t>(j*nx_+i)*ng;
            geometry_[q+3]=geometry_[static_cast<std::size_t>(i)*ng+3];
            geometry_[q+5]=geometry_[static_cast<std::size_t>(i)*ng+5];
            geometry_[q+4]=geometry_[static_cast<std::size_t>((ny_-1)*nx_+i)*ng+4];
            if(std::abs(geometry_[q+4]-3.035)>1e-9) throw std::invalid_argument("hill top wall must be y=3.035");
        }
    }
    if(std::abs(crest_area-area)>1e-9) throw std::invalid_argument("hill crest area differs");
}

void PeriodicHillRuntime::bind_force(std::function<void(Real)> setter) {
    setter_=std::move(setter); setter_(force_);
}
void PeriodicHillRuntime::sample_flow() {
    std::fill(instant_.begin(),instant_.end(),0);
    const auto& c=context_;
    bool valid=true;
    for(const auto& b:c.local_blocks.blocks()) {
        const auto& l=leaf_for(c,b.id()); const auto e=b.cell_extent();
        for(int k=0;k<e.nk;++k) for(int j=0;j<e.nj;++j) for(int i=0;i<e.ni;++i) {
            const Real dz=b.coordinates.z(i,j,k+1)-b.coordinates.z(i,j,k);
            const auto q=static_cast<std::size_t>((j+l.cells.begin.j)*nx_+i+l.cells.begin.i)*nm;
            const auto p=temperature_primitive_from_conservative(load_conservative(b.flow.conservative,{i,j,k}),
                c.quantities.gas,c.quantities.reference,c.quantities.floors,3);
            const auto pp=pressure_primitive(p,c.quantities.gas,c.quantities.reference,c.quantities.floors,3);
            const Real r=p[0],u=p[1],v=p[2],w=p[3];
            const Real a[nm]={r,u,v,w,pp[4],p[4],u*u,v*v,w*w,u*v,u*w,v*w,r*u};
            for(int m=0;m<nm;++m) { valid=valid && std::isfinite(a[m]); instant_[q+m]+=dz*a[m]; }
        }
    }
    if(!c.mpi.all_true(valid)) throw std::runtime_error("hill: nonfinite sample");
    c.mpi.sum_reals(instant_);
    for(std::size_t n=0;n<instant_.size()/nm;++n)
        for(int m=0;m<nm;++m) instant_[n*nm+m]/=geometry_[n*ng+2];
    mass_flux_=volume_flux_=0;
    // Periodic seam interpolation, not a sample at the first positive-x center.
    const Real x0=geometry_[0], x1=geometry_[static_cast<std::size_t>(nx_-1)*ng];
    const Real a=(9-x1)/(9-x1+x0);
    for(int j=0;j<ny_;++j) {
        const auto l=static_cast<std::size_t>(j*nx_),r=l+nx_-1;
        const Real A=geometry_[l*ng+6];
        mass_flux_+=A*(a*instant_[l*nm+12]+(1-a)*instant_[r*nm+12]);
        volume_flux_+=A*(a*instant_[l*nm+1]+(1-a)*instant_[r*nm+1]);
    }
    const Real mu=c.quantities.transport.viscosity(1)/c.quantities.reference.reynolds();
    for(int i=0;i<nx_;++i) for(int side=0;side<2;++side) {
        const int j1=side ? ny_-1 : 0, j2=side ? ny_-2 : 1;
        const auto n1=static_cast<std::size_t>(j1*nx_+i), n2=static_cast<std::size_t>(j2*nx_+i);
        const Real slope=side ? 0 : geometry_[n1*ng+5];
        const Real norm=std::sqrt(1+slope*slope);
        const Real wall_y=geometry_[n1*ng+(side ? 4 : 3)];
        const Real d1=std::abs(geometry_[n1*ng+1]-wall_y),d2=std::abs(geometry_[n2*ng+1]-wall_y);
        const Real u1=(instant_[n1*nm+1]+slope*instant_[n1*nm+2])/norm;
        const Real u2=(instant_[n2*nm+1]+slope*instant_[n2*nm+2])/norm;
        wall_[static_cast<std::size_t>(i)*4+side*2]=mu*norm*hill_wall_derivative(u1,u2,d1,d2);
        // Linear pressure extrapolation to the wall.
        wall_[static_cast<std::size_t>(i)*4+side*2+1]
            =(d2*instant_[n1*nm+4]-d1*instant_[n2*nm+4])/(d2-d1);
    }
}
void PeriodicHillRuntime::on_initial(const SimulationState& s) {
    sample_flow();
    if(!restored_) { previous_error_=config_.periodic_hill.bulk_velocity-mass_flux_/(config_.periodic_hill.density*area); last_step_=s.step; last_time_=s.time; }
    if(last_step_!=s.step || std::abs(last_time_-s.time)>1e-12)
        throw std::runtime_error("hill checkpoint step/time mismatch");
}
void PeriodicHillRuntime::on_step(const SimulationState& s,bool) {
    if(!s.residuals.finite || s.step<=last_step_ || s.time<=last_time_) return;
    sample_flow();
    const auto& h=config_.periodic_hill;
    const Real dt=s.time-last_time_;
    const Real error=h.bulk_velocity-mass_flux_/(h.density*area);
    force_=hill_feedback_force(force_,error,previous_error_,dt,h.controller_time,h.density,h.force_limit);
    previous_error_=error;
    if(setter_) setter_(force_);
    // Right-endpoint quadrature of each accepted physical interval, clipped to the window.
    const Real weight=std::max(Real(0),std::min(s.time,h.statistics_end)-std::max(last_time_,h.statistics_start));
    if(weight>0) {
        for(std::size_t q=0;q<integral_.size();++q) integral_[q]+=weight*instant_[q];
        for(std::size_t q=0;q<wall_.size();++q) wall_integral_[q]+=weight*wall_[q];
        weight_+=weight; ++samples_;
    }
    last_step_=s.step; last_time_=s.time;
    if(s.step%h.write_every_steps==0) { write_history(s); write_averages(); }
}
void PeriodicHillRuntime::on_final(const SimulationState& s) {
    // A failed initial residual can reach on_final before output-directory setup.
    // Keep the last successful samples; the output manager writes the failure manifest.
    if(s.stop_reason==StopReason::NumericalFailure) return;
    if(s.step!=last_history_step_) write_history(s);
    write_averages();
}
std::vector<std::string> PeriodicHillRuntime::output_paths() const {
    std::vector<std::string> out;
    for(const auto* suffix:{"hill_history.csv","hill_mean_xy.csv","hill_profiles.csv","hill_wall.csv"})
        out.push_back(config_.output.directory+"/"+config_.case_name+"_"+suffix);
    return out;
}
void PeriodicHillRuntime::write_history(const SimulationState& s) {
    root_io(context_.mpi,[&] {
        const auto path=output_paths()[0]; struct stat info {};
        const bool exists=stat(path.c_str(), &info)==0;
        auto f=output_file(path,exists);
        if(!exists) f<<"step,time,force_next,mass_flux,volume_flux,Ub_mass,Ub_volume,Re_mass,statistics_weight,samples\n";
        const Real ub=mass_flux_/(config_.periodic_hill.density*area);
        f<<s.step<<','<<s.time<<','<<force_<<','<<mass_flux_<<','<<volume_flux_<<','<<ub<<','<<volume_flux_/area<<','
         <<config_.periodic_hill.density*ub*context_.quantities.reference.reynolds()/context_.quantities.transport.viscosity(1)
         <<','<<weight_<<','<<samples_<<'\n';
    });
    last_history_step_=s.step;
}
void PeriodicHillRuntime::write_averages() const {
    root_io(context_.mpi,[&] {
        const auto paths=output_paths(); auto xy=output_file(paths[1]),pr=output_file(paths[2]),wall=output_file(paths[3]);
        const char* columns="x_h,y_h,rho,u_Ub,v_Ub,w_Ub,p,T,uu_Ub2,vv_Ub2,ww_Ub2,uv_Ub2,uw_Ub2,vw_Ub2\n";
        xy<<columns; pr<<columns; wall<<"x_h,lower_y_h,tau_lower,Cf_lower,p_lower,tau_upper,Cf_upper,p_upper\n";
        if(weight_<=0) return; // Empty table explicitly means no accepted statistics yet.
        const Real U=config_.periodic_hill.bulk_velocity;
        auto mean=[&](std::size_t q) {
            std::vector<Real> m(nm);
            for(int a=0;a<nm;++a) m[a]=integral_[q*nm+a]/weight_;
            const int c1[6]={1,2,3,1,1,2},c2[6]={1,2,3,2,3,3};
            for(int a=0;a<6;++a) {
                const Real raw=m[6+a];
                m[6+a]-=m[c1[a]]*m[c2[a]];
                if(a<3 && m[6+a]<0 && m[6+a]>-1e-12*std::max(Real(1),raw)) m[6+a]=0;
            }
            return m;
        };
        auto row=[&](std::ostream& f,Real x,Real y,const std::vector<Real>& m) {
            f<<x<<','<<y<<','<<m[0]<<','<<m[1]/U<<','<<m[2]/U<<','<<m[3]/U<<','<<m[4]<<','<<m[5];
            for(int a=0;a<6;++a) f<<','<<m[6+a]/(U*U);
            f<<'\n';
        };
        for(int j=0;j<ny_;++j) for(int i=0;i<nx_;++i) {
            const auto q=static_cast<std::size_t>(j*nx_+i); const auto m=mean(q);
            row(xy,geometry_[q*ng],geometry_[q*ng+1],m);
        }
        // Form covariance at each cell before spatial interpolation: a laminar
        // streamwise gradient must not become a spurious Reynolds stress.
        for(Real x:{.05,2.,4.,8.}) {
            int hi=0; while(hi<nx_ && geometry_[static_cast<std::size_t>(hi)*ng]<x) ++hi;
            const int lo=(hi+nx_-1)%nx_; hi%=nx_;
            Real xl=geometry_[static_cast<std::size_t>(lo)*ng],xr=geometry_[static_cast<std::size_t>(hi)*ng];
            if(x<xl) xl-=9;
            if(x>xr) xr+=9;
            const Real a=(x-xl)/(xr-xl);
            for(int j=0;j<ny_;++j) {
                const auto l=static_cast<std::size_t>(j*nx_+lo),r=static_cast<std::size_t>(j*nx_+hi);
                std::vector<Real> m(nm); const auto ml=mean(l),mr=mean(r);
                for(int v=0;v<nm;++v) m[v]=(1-a)*ml[v]+a*mr[v];
                row(pr,x,(1-a)*geometry_[l*ng+1]+a*geometry_[r*ng+1],m);
            }
        }
        const Real dynamic=.5*config_.periodic_hill.density*U*U;
        for(int i=0;i<nx_;++i) {
            const auto q=static_cast<std::size_t>(i);
            wall<<geometry_[q*ng]<<','<<geometry_[q*ng+3];
            for(int side=0;side<2;++side) {
                const Real tau=wall_integral_[q*4+side*2]/weight_;
                wall<<','<<tau<<','<<tau/dynamic<<','<<wall_integral_[q*4+side*2+1]/weight_;
            }
            wall<<'\n';
        }
    });
}
std::string PeriodicHillRuntime::serialize() const {
    std::ostringstream s; s<<std::setprecision(17)<<config_.periodic_hill.signature()<<'\n';
    s<<nx_<<' '<<ny_<<' '<<force_<<' '<<previous_error_<<' '<<weight_<<' '<<samples_<<' '<<last_step_<<' '<<last_time_<<'\n';
    for(Real v:integral_) s<<v<<' ';
    s<<'\n';
    for(Real v:wall_integral_) s<<v<<' ';
    s<<'\n';
    return s.str();
}
void PeriodicHillRuntime::restore(const std::string& text) {
    std::istringstream s(text); std::string signature; std::getline(s,signature);
    if(signature!=config_.periodic_hill.signature()) throw std::runtime_error("hill checkpoint signature differs");
    int nx=0,ny=0; Real force=0,error=0,weight=0,time=0; std::size_t samples=0,step=0;
    s>>nx>>ny>>force>>error>>weight>>samples>>step>>time;
    if(!s || nx!=nx_ || ny!=ny_ || !std::isfinite(force) || !std::isfinite(error)
        || !std::isfinite(weight) || weight<0 || !std::isfinite(time) || time<0
        || std::abs(force)>config_.periodic_hill.force_limit)
        throw std::runtime_error("invalid hill checkpoint header");
    auto integral=integral_,wall=wall_integral_;
    for(auto* values:{&integral,&wall}) for(auto& v:*values) {
        s>>v; if(!s || !std::isfinite(v)) throw std::runtime_error("invalid hill checkpoint moments");
    }
    std::string extra; if(s>>extra) throw std::runtime_error("trailing hill checkpoint data");
    force_=force; previous_error_=error; weight_=weight; samples_=samples; last_step_=step; last_time_=time;
    integral_=std::move(integral); wall_integral_=std::move(wall); restored_=true;
}
} // namespace wcns
