#include <wcns/runtime/hit.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>


namespace wcns {
namespace {
constexpr Real pi=3.141592653589793238462643383279502884;
const std::vector<std::string> names={
 "rho_mean","u_mean","v_mean","w_mean","p_mean","T_mean",
 "uu","vv","ww","uv","uw","vw","favre_u","favre_v","favre_w",
 "favre_uu","favre_vv","favre_ww","favre_uv","favre_uw","favre_vw",
 "rho_rms","p_rms","T_rms","K","K_mass","urms","epsilon_viscous",
 "enstrophy","divergence_rms","pressure_dilatation","helicity","L_integral",
 "lambda_isotropic","Re_lambda_isotropic","eta","kmax_eta","Mt",
 "K_solenoidal","K_dilatational","K_forced_band","forcing_power","cooling",
 "forcing_alpha","du_dx_skewness","du_dx_flatness","dv_dy_skewness","dv_dy_flatness",
 "dw_dz_skewness","dw_dz_flatness","parseval_error","mass","total_energy"};
constexpr std::size_t Kmindex=25, Epsindex=27, Pdindex=30, Pinindex=41;
bool exists(const std::string& p) { return std::ifstream(p).good(); }
Real sq(Real x) {return x*x;}
std::uint64_t mix(std::uint64_t x) {
    x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
}
Real random_signed(std::uint64_t x) {return 2*static_cast<Real>(mix(x)>>11)*0x1.0p-53-1;}
std::ofstream output(const std::string& path,bool append=false) {
    std::ofstream f(path,append?std::ios::app:std::ios::out);
    f.exceptions(std::ios::badbit|std::ios::failbit);f<<std::setprecision(17);return f;
}
Real interpolate(const std::vector<std::array<Real,2>>& table, Real k) {
    if(k<=0) return 0;
    if(k<table.front()[0]) return table.front()[1]*std::pow(k/table.front()[0],4);
    if(k>table.back()[0]) return 0;
    auto it=std::lower_bound(table.begin(),table.end(),k,[](const auto& a,Real x){return a[0]<x;});
    if(it==table.begin()) return (*it)[1];
    const auto& a=*(it-1);const auto& b=*it;
    const Real t=std::log(k/a[0])/std::log(b[0]/a[0]);
    if(a[1]==0 || b[1]==0) return a[1]+t*(b[1]-a[1]);
    return std::exp(std::log(a[1])+t*std::log(b[1]/a[1]));
}
Real kang_spectrum(const std::vector<std::array<Real,2>>& table, Real k) {
    // energy.m extrapolates the first/last log-log segment; the legacy
    // initializer's k^4 low-k continuation is deliberately kept separate.
    if(k>=table.front()[0] && k<=table.back()[0])return interpolate(table,k);
    const std::size_t i=k<table.front()[0]?0:table.size()-2;
    const auto& a=table[i];const auto& b=table[i+1];
    if(a[1]==0 || b[1]==0)return 0;
    return a[1]*std::pow(k/a[0],std::log(b[1]/a[1])/std::log(b[0]/a[0]));
}
}
void HitConfig::validate() const {
    if(!enabled()) return;
    if((initialization!="shell_spectrum" && initialization!="analytic_random_phase")
       || !std::isfinite(spectrum_amplitude) || spectrum_amplitude<=0
       || !std::isfinite(preparation_time) || preparation_time<0 || (preparation_time>0 && preparation_time<1e-12) || !preparation_max_steps)
        throw CaseConfigurationError("invalid HIT initialization, spectrum amplitude or preparation settings");
    if(initialization=="analytic_random_phase" && (!spectrum_file.empty() || (type=="forced" && forcing=="jhtdb_shells")))
        throw CaseConfigurationError("analytic_random_phase requires no spectrum_file and no JHTDB shell override");
    if(preparation_time>0 && (type!="decay" || initialization!="shell_spectrum" || spectrum_file.empty()))
        throw CaseConfigurationError("HIT preparation requires decay, shell_spectrum and a target spectrum_file");
    if((type!="decay" && type!="forced") || n<8 || n>512 || (n&(n-1))
       || !std::isfinite(length) || length<=0 || !std::isfinite(initial_energy) || initial_energy<=0
       || !std::isfinite(peak_wave) || peak_wave<=0 || !std::isfinite(cutoff) || cutoff<0 || cutoff>=n/2
       || (forcing!="constant_power" && forcing!="constant_band_energy" && forcing!="jhtdb_shells")
       || !std::isfinite(forcing_power) || forcing_power<=0 || !std::isfinite(forcing_kmax)
       || forcing_kmax<1 || forcing_kmax>=n/2 || !sample_every_steps || !write_every_samples
       || !std::isfinite(statistics_start) || !std::isfinite(statistics_end)
       || statistics_start<0 || statistics_end<=statistics_start)
        throw CaseConfigurationError("invalid HIT type, cube, spectrum, forcing or statistics settings");
    if(forcing=="jhtdb_shells" && type=="forced" && std::abs(length-2*pi)>1e-12)
        throw CaseConfigurationError("JHTDB shells require box length 2*pi");
    if(forcing=="jhtdb_shells" && type=="forced" && (cutoff>0?cutoff:Real(n)/3)<2.5)
        throw CaseConfigurationError("JHTDB initialization must include both complete forced shells (cutoff >= 2.5)");
    Real previous=-1;
    for(Real t:sample_times) {
        if(!std::isfinite(t) || t<0 || t<=previous)
            throw CaseConfigurationError("hit.sample_times must be nonnegative and strictly increasing");
        previous=t;
    }
}
std::string HitConfig::signature() const {
    std::ostringstream s;s<<std::setprecision(17)<<"hit_v1;"<<type<<';'<<n<<';'<<length<<';'<<seed
      <<';'<<spectrum_file<<';'<<initial_energy<<';'<<peak_wave<<';'<<cutoff<<';'<<forcing<<';'<<forcing_power
      <<';'<<forcing_kmax<<';'<<remove_mean_acceleration<<';'<<thermostat<<';'<<statistics_start
      <<';'<<statistics_end<<';'<<sample_every_steps;
    for(Real t:sample_times) s<<';'<<t;
    // Preserve strict restart compatibility with v2.6 for unchanged legacy cases.
    if(initialization!="shell_spectrum" || preparation_time>0)
        s<<";init_v2;method="<<initialization<<";A="<<spectrum_amplitude<<";preparation="<<preparation_time;
    return s.str();
}
HitRuntime::HitRuntime(const CaseConfig& config,const StatisticContext& context,
                       LocalBlockSet& blocks,std::string spectrum_path)
    :config_(config),context_(context),blocks_(blocks),fft_(context.mpi,config.hit.n)
{
    const auto& mpi=context.mpi;const int n=fft_.n();const Real dx=config.hit.length/n;
    bool valid=true;
    // Strict geometry contract: uniform, axis-aligned periodic cube [0,L]^3.
    for(auto& b:blocks_.blocks()) {
        valid=valid && b.cell_dimension()==3;
        for(const auto& patch:b.boundaries) valid=valid && patch.type==BoundaryType::Periodic;
        const auto& m=context.metrics.at(b.id());const auto e=b.cell_extent();
        for(int k=0;k<e.nk;++k) for(int j=0;j<e.nj;++j) for(int i=0;i<e.ni;++i) {
            std::array<Real,3> x{{m.cell_coordinates().x(i,j,k),m.cell_coordinates().y(i,j,k),m.cell_coordinates().z(i,j,k)}};
            std::array<int,3> a{};
            for(int d=0;d<3;++d) {a[d]=static_cast<int>(std::llround(x[d]/dx-.5));
                valid=valid && a[d]>=0 && a[d]<n && std::abs(x[d]-(a[d]+.5)*dx)<1e-8*dx;}
            valid=valid && std::abs(m.jacobian()(i,j,k)-dx*dx*dx)<1e-7*dx*dx*dx;
            cells_.push_back({&b,{i,j,k},(static_cast<std::size_t>(std::max(0,a[0]))*n+std::max(0,a[1]))*n+std::max(0,a[2])});
        }
    }
    if(!mpi.all_true(valid)) throw std::runtime_error("HIT requires uniform Cartesian cube cells and no physical boundary patches");
    send_counts_.resize(mpi.size());std::vector<Real> indices;
    for(int r=0;r<mpi.size();++r) for(std::size_t c=0;c<cells_.size();++c)
        if(fft_.owner(static_cast<int>(cells_[c].global/(n*n)))==r) {
            ++send_counts_[r];send_order_.push_back(c);indices.push_back(static_cast<Real>(cells_[c].global));
        }
    auto received=exchange_peer_reals(mpi,indices,send_counts_,&receive_counts_);
    std::vector<int> seen(fft_.size());
    for(Real q:received) {
        const auto index=static_cast<std::size_t>(q)-static_cast<std::size_t>(fft_.begin())*n*n;
        if(index>=seen.size()) {valid=false;received_indices_.push_back(0);}
        else {valid=valid && ++seen[index]==1;received_indices_.push_back(index);}
    }
    valid=valid && std::all_of(seen.begin(),seen.end(),[](int v){return v==1;});
    if(!mpi.all_true(valid)) throw std::runtime_error("HIT cube has duplicate or missing cells");
    std::string content,error;
    if(mpi.rank()==0 && !spectrum_path.empty()) {
        try {std::ifstream f(spectrum_path);if(!f) throw std::runtime_error("cannot read HIT spectrum "+spectrum_path);
            content.assign(std::istreambuf_iterator<char>(f),{});}catch(const std::exception& e){error=e.what();}
    }
    error=mpi.broadcast_string(error);if(!error.empty()) throw std::runtime_error(error);
    content=mpi.broadcast_string(content);
    std::uint64_t hash=14695981039346656037ULL;
    for(unsigned char ch:content) {hash^=ch;hash*=1099511628211ULL;}
    spectrum_identity_=std::to_string(hash);
    if(!spectrum_path.empty()) {
        std::istringstream in(content);std::string line;Real last=0;
        while(std::getline(in,line)) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream row(line);Real k,e;std::string extra;
            if(!(row>>k>>e) || (row>>extra) || !std::isfinite(k) || !std::isfinite(e) || k<=last || e<0)
                throw std::runtime_error("HIT spectrum requires increasing positive k and nonnegative E(k), two whitespace columns");
            target_spectrum_.push_back({k,e});last=k;
        }
        if(target_spectrum_.size()<2) throw std::runtime_error("HIT spectrum is empty");
    }
    instant_.resize(names.size());previous_=integral_=instant_;
    spectrum_integral_.resize(static_cast<std::size_t>(std::ceil(std::sqrt(3.)*n/2))+1);
    preparation_complete_=config.hit.preparation_time==0;
}
std::vector<Real> HitRuntime::to_slab(const std::vector<Real>& v,int nc) const {
    if(v.size()!=cells_.size()*nc) throw std::invalid_argument("HIT cell payload");
    std::vector<Real> send;send.reserve(v.size());
    for(auto q:send_order_) for(int c=0;c<nc;++c) send.push_back(v[q*nc+c]);
    auto counts=send_counts_;for(auto& c:counts)c*=nc;
    auto recv=exchange_peer_reals(context_.mpi,send,counts);
    std::vector<Real> out(fft_.size()*nc);
    for(std::size_t q=0;q<received_indices_.size();++q) for(int c=0;c<nc;++c)
        out[received_indices_[q]*nc+c]=recv[q*nc+c];
    return out;
}
std::vector<Real> HitRuntime::from_slab(const std::vector<Real>& v,int nc) const {
    if(v.size()!=fft_.size()*nc) throw std::invalid_argument("HIT slab payload");
    std::vector<Real> send;send.reserve(v.size());
    for(auto q:received_indices_) for(int c=0;c<nc;++c) send.push_back(v[q*nc+c]);
    auto counts=receive_counts_;for(auto& c:counts)c*=nc;
    auto recv=exchange_peer_reals(context_.mpi,send,counts);
    std::vector<Real> out(cells_.size()*nc);
    for(std::size_t q=0;q<send_order_.size();++q) for(int c=0;c<nc;++c)
        out[send_order_[q]*nc+c]=recv[q*nc+c];
    return out;
}
std::vector<Real> HitRuntime::velocities(bool residual) const {
    std::vector<Real> v(cells_.size()*3);
    for(std::size_t q=0;q<cells_.size();++q) {const auto& c=cells_[q];const auto a=c.index;
        const auto u=load_conservative(c.block->flow.conservative,a);
        for(int d=0;d<3;++d) v[q*3+d]=residual?
            (c.block->flow.residual(a.i,a.j,a.k,d+1)-u[d+1]/u[0]*c.block->flow.residual(a.i,a.j,a.k,0))/u[0]:u[d+1]/u[0];
    }return to_slab(v,3);
}
std::array<FourierField,3> HitRuntime::velocity_spectrum() const {
    auto v=velocities();std::array<FourierField,3> u;
    for(int d=0;d<3;++d) {u[d].resize(fft_.size());for(std::size_t q=0;q<fft_.size();++q)u[d][q]=v[3*q+d];fft_.forward(u[d]);}
    return u;
}
void HitRuntime::initialize() {
    const auto& h=config_.hit;const int n=fft_.n();const Real dk=2*pi/h.length;
    const Real cutoff=h.cutoff>0?h.cutoff:Real(n)/3;
    const bool analytic=h.initialization=="analytic_random_phase";
    std::array<FourierField,3> u;for(auto& a:u)a.resize(fft_.size());
    std::vector<Real> shell(spectrum_integral_.size());
    for(std::size_t q=0;q<fft_.size();++q) {
        const auto k=fft_.wave(q);const Real k2=sq(k[0])+sq(k[1])+sq(k[2]);
        if(!k2 || std::sqrt(k2)>cutoff) continue;
        const int first=k[0]!=0?k[0]:(k[1]!=0?k[1]:k[2]);
        const int sign=first>0?1:-1;
        const auto canonical=(static_cast<std::uint64_t>(sign*k[0]+512)*1025+(sign*k[1]+512))*1025+(sign*k[2]+512);
        FourierValue dot=0;
        for(int d=0;d<3;++d) {
            if(analytic) {
                // Random phases, transverse projection, then ONE global energy
                // normalization. Canonical +/- modes ensure a real inverse FFT.
                const Real phase=pi*(random_signed(canonical*7+d+h.seed)+1);
                const Real physical_k2=k2*dk*dk;
                const Real amplitude=std::sqrt((2./3)*h.spectrum_amplitude*physical_k2
                    *std::exp(-2*physical_k2/sq(h.peak_wave))*dk*dk*dk/(4*pi));
                u[d][q]=amplitude*FourierValue(std::sin(phase),sign*std::cos(phase));
            } else {
                u[d][q]={random_signed(canonical*7+2*d+h.seed),sign*random_signed(canonical*7+2*d+1+h.seed)};
            }
            dot+=Real(k[d])*u[d][q];
        }
        Real e=0;for(int d=0;d<3;++d) {u[d][q]-=Real(k[d])/k2*dot;e+=.5*std::norm(u[d][q]);}
        shell[static_cast<std::size_t>(std::floor(std::sqrt(k2)+.5))]+=e;
    }
    context_.mpi.sum_reals(shell);
    if(analytic) {
        const Real energy=std::accumulate(shell.begin(),shell.end(),Real{0});
        const Real desired=3*h.spectrum_amplitude/64*std::sqrt(2*pi)*std::pow(h.peak_wave,5);
        if(!(energy>0) || !std::isfinite(desired) || !(desired>0))
            throw std::runtime_error("HIT analytic spectrum has no finite resolved energy");
        const Real scale=std::sqrt(desired/energy);
        for(auto& a:u)for(auto& v:a)v*=scale;
        install_velocity(std::move(u),true);
        return;
    }
    std::vector<Real> target(shell.size());Real total=0;
    for(std::size_t s=1;s<shell.size();++s) if(shell[s]>0) {
        const Real k=static_cast<Real>(s)*dk;
        target[s]=target_spectrum_.empty()?std::pow(k,4)*std::exp(-2*sq(k/h.peak_wave))*dk:interpolate(target_spectrum_,k)*dk;
        total+=target[s];
    }
    if(!(total>0)) throw std::runtime_error("HIT initial spectrum has no resolved energy");
    // Experimental spectra are NOT renormalized to put unresolved energy on coarse modes.
    if(target_spectrum_.empty()) for(auto& e:target)e*=h.initial_energy/total;
    if(h.type=="forced" && h.forcing=="jhtdb_shells") {target[1]=.30;target[2]=.13;}
    for(std::size_t q=0;q<fft_.size();++q) {
        auto k=fft_.wave(q);auto s=static_cast<std::size_t>(std::floor(std::sqrt(sq(k[0])+sq(k[1])+sq(k[2]))+.5));
        const Real scale=shell[s]>0?std::sqrt(target[s]/shell[s]):0;
        for(auto& a:u)a[q]*=scale;
    }
    install_velocity(std::move(u),true);
}
void HitRuntime::install_velocity(std::array<FourierField,3> u,bool fresh) {
    for(auto& a:u)fft_.inverse(a);
    std::vector<Real> v(fft_.size()*3);Real imaginary=0;
    for(std::size_t q=0;q<fft_.size();++q) for(int d=0;d<3;++d) {v[q*3+d]=u[d][q].real();imaginary=std::max(imaginary,std::abs(u[d][q].imag()));}
    if(context_.mpi.max(imaginary)>1e-10) throw std::runtime_error("HIT initial Hermitian symmetry failed");
    auto local=from_slab(v,3);
    const auto& gas=context_.quantities.gas;const auto& ref=context_.quantities.reference;
    const Real internal=1/(gas.gamma()*(gas.gamma()-1)*sq(ref.mach()));
    for(std::size_t q=0;q<cells_.size();++q) {auto& c=cells_[q];const auto a=c.index;
        const auto old=load_conservative(c.block->flow.conservative,a);
        const Real rho=fresh?1:old[0];
        Real energy=fresh?internal:old[4]-.5*(sq(old[1])+sq(old[2])+sq(old[3]))/rho;
        c.block->flow.conservative(a.i,a.j,a.k,0)=rho;
        for(int d=0;d<3;++d) {const Real v0=local[q*3+d];c.block->flow.conservative(a.i,a.j,a.k,d+1)=rho*v0;energy+=.5*rho*v0*v0;}
        c.block->flow.conservative(a.i,a.j,a.k,4)=energy;
    }
}
void HitRuntime::finish_preparation(const SimulationState& state) {
    const auto& h=config_.hit;
    if(!preparing() || state.step==0 || state.stop_reason!=StopReason::PhysicalTimeReached
       || std::abs(state.time-h.preparation_time)>1e-12*std::max(Real{1},h.preparation_time))
        throw std::runtime_error("HIT spectrum rematch requires a completed preparation stage");
    auto u=velocity_spectrum();
    const Real dk=2*pi/h.length,cutoff=h.cutoff>0?h.cutoff:Real(fft_.n())/3;
    std::vector<Real> audit(spectrum_integral_.size()*3);bool valid=true;
    for(std::size_t q=0;q<fft_.size();++q) {
        const auto m=fft_.wave(q);const Real m2=sq(m[0])+sq(m[1])+sq(m[2]),km=std::sqrt(m2);
        const auto shell=static_cast<std::size_t>(std::floor(km+.5));
        for(auto& a:u)audit[3*shell]+=.5*std::norm(a[q]);
        if(!m2 || km>cutoff) {for(auto& a:u)a[q]=0;continue;}
        // The author's incompressible input is already transverse. Our weakly
        // compressible preparation needs projection before Kang rescaling.
        FourierValue dot=0;for(int d=0;d<3;++d)dot+=Real(m[d])*u[d][q];
        Real magnitude2=0;
        for(int d=0;d<3;++d){u[d][q]-=Real(m[d])/m2*dot;magnitude2+=std::norm(u[d][q]);}
        const Real desired2=2*kang_spectrum(target_spectrum_,km*dk)*dk/(4*pi*m2);
        valid=valid && std::isfinite(desired2) && (desired2==0 || magnitude2>1e-28*desired2);
        const Real factor=magnitude2>0?std::sqrt(desired2/magnitude2):0;
        for(auto& a:u){a[q]*=factor;audit[3*shell+1]+=.5*std::norm(a[q]);}
        audit[3*shell+2]+=.5*desired2;
    }
    if(!context_.mpi.all_true(valid))throw std::runtime_error("HIT rematch has a zero transverse mode with nonzero target energy");
    context_.mpi.sum_reals(audit);
    root_io([&]{auto f=output(config_.output.directory+"/preparation/"+config_.case_name+"_rematch.csv");
        f<<"shell,k,E_before,E_after,E_target_discrete\n";
        for(std::size_t s=0;s<audit.size()/3;++s)f<<s<<','<<s*dk<<','<<audit[3*s]/dk<<','<<audit[3*s+1]/dk<<','<<audit[3*s+2]/dk<<'\n';});
    install_velocity(std::move(u),false); // preserve evolved density and internal energy
    preparation_complete_=true;preparation_steps_=state.step;
    std::fill(instant_.begin(),instant_.end(),0);previous_=integral_=instant_;
    std::fill(spectrum_integral_.begin(),spectrum_integral_.end(),0);
    weight_=last_sample_time_=last_step_time_=force_power_=cooling_=force_alpha_=0;
    samples_=last_sample_step_=last_step_=0;restored_=sampled_=false;
}
void HitRuntime::add_stage_source(Real) {
    if(config_.hit.forcing!="jhtdb_shells")force_power_=0;
    cooling_=force_alpha_=0;
    if(config_.hit.type!="forced") return;
    if(config_.hit.forcing=="jhtdb_shells") {
        Real rate=0;
        for(const auto& c:cells_) {const auto b=c.index;const auto u0=load_conservative(c.block->flow.conservative,b);
            Real ke=0,ri=c.block->flow.residual(b.i,b.j,b.k,4);
            for(int d=0;d<3;++d) {const Real v=u0[d+1]/u0[0];ke+=.5*v*v;ri-=v*c.block->flow.residual(b.i,b.j,b.k,d+1);}
            rate+=ri+ke*c.block->flow.residual(b.i,b.j,b.k,0);
        }
        const Real count=static_cast<Real>(fft_.n())*fft_.n()*fft_.n();
        cooling_=config_.hit.thermostat?context_.mpi.sum(rate)/count:0;
        for(const auto& c:cells_){const auto b=c.index;c.block->flow.residual(b.i,b.j,b.k,4)-=cooling_;}
        return;
    }
    auto u=velocity_spectrum();auto a=u;Real band=0,rate=0;
    std::array<FourierField,3> du;
    const bool constant_band=config_.hit.forcing=="constant_band_energy";
    if(constant_band) {auto v=velocities(true);for(int d=0;d<3;++d){du[d].resize(fft_.size());
        for(std::size_t q=0;q<fft_.size();++q)du[d][q]=v[3*q+d];
        fft_.forward(du[d]);}}
    for(std::size_t q=0;q<fft_.size();++q) {
        const auto k=fft_.wave(q);const Real k2=sq(k[0])+sq(k[1])+sq(k[2]);
        FourierValue dot=0;for(int d=0;d<3;++d)dot+=Real(k[d])*u[d][q];
        for(int d=0;d<3;++d) {
            a[d][q]=(k2>0 && k2<=sq(config_.hit.forcing_kmax))?u[d][q]-Real(k[d])*dot/k2:FourierValue{};
            band+=std::norm(a[d][q]);
            if(constant_band)rate+=(std::conj(a[d][q])*du[d][q]).real();
        }
    }
    band=context_.mpi.sum(band);rate=context_.mpi.sum(rate);
    if(!(band>1e-20)) throw std::runtime_error("HIT forced band has zero energy");
    for(auto& f:a)fft_.inverse(f);
    std::vector<Real> slab(fft_.size()*3);
    for(std::size_t q=0;q<fft_.size();++q)for(int d=0;d<3;++d)slab[q*3+d]=a[d][q].real();
    auto acc=from_slab(slab,3);std::vector<Real> sums(4);
    for(std::size_t q=0;q<cells_.size();++q) {const auto& c=cells_[q];const auto u0=load_conservative(c.block->flow.conservative,c.index);
        sums[0]+=u0[0];for(int d=0;d<3;++d)sums[d+1]+=u0[0]*acc[q*3+d];}
    context_.mpi.sum_reals(sums);
    Real work=0,internal_rate=0;
    for(std::size_t q=0;q<cells_.size();++q) {const auto& c=cells_[q];const auto u0=load_conservative(c.block->flow.conservative,c.index);
        const auto b=c.index;Real ke=0,ri=c.block->flow.residual(b.i,b.j,b.k,4);
        for(int d=0;d<3;++d) {
            if(config_.hit.remove_mean_acceleration)acc[q*3+d]-=sums[d+1]/sums[0];
            work+=u0[d+1]*acc[q*3+d];ke+=.5*sq(u0[d+1]/u0[0]);
            ri-=u0[d+1]/u0[0]*c.block->flow.residual(b.i,b.j,b.k,d+1);
        }internal_rate+=ri+ke*c.block->flow.residual(b.i,b.j,b.k,0);
    }
    work=context_.mpi.sum(work);
    if(!constant_band && !(work>1e-20))throw std::runtime_error("HIT forcing work is nonpositive");
    force_alpha_=constant_band?-rate/band:config_.hit.forcing_power*sums[0]/work;
    force_power_=force_alpha_*work/sums[0];
    const Real count=static_cast<Real>(fft_.n())*fft_.n()*fft_.n();
    cooling_=config_.hit.thermostat?context_.mpi.sum(internal_rate)/count:0;
    for(std::size_t q=0;q<cells_.size();++q) {const auto& c=cells_[q];const auto b=c.index;
        const auto u0=load_conservative(c.block->flow.conservative,b);Real power=0;
        for(int d=0;d<3;++d) {const Real f=force_alpha_*acc[q*3+d];
            c.block->flow.residual(b.i,b.j,b.k,d+1)+=u0[0]*f;power+=u0[d+1]*f;}
        c.block->flow.residual(b.i,b.j,b.k,4)+=power-cooling_;
    }
}

void HitRuntime::accepted_step_transform(Real dt) {
    if(config_.hit.type!="forced" || config_.hit.forcing!="jhtdb_shells")return;
    if(!(dt>0))throw std::invalid_argument("HIT rescale time step");
    auto u=velocity_spectrum();std::vector<Real> energy(2);
    for(std::size_t q=0;q<fft_.size();++q){auto k=fft_.wave(q);Real k2=sq(k[0])+sq(k[1])+sq(k[2]);
        int shell=static_cast<int>(std::floor(std::sqrt(k2)+.5));if(shell<1 || shell>2)continue;
        FourierValue dot=0;for(int d=0;d<3;++d)dot+=Real(k[d])*u[d][q];
        for(int d=0;d<3;++d)energy[shell-1]+=.5*std::norm(u[d][q]-Real(k[d])*dot/k2);
    }
    context_.mpi.sum_reals(energy);
    if(energy[0]<=1e-20 || energy[1]<=1e-20)throw std::runtime_error("JHTDB shell energy vanished");
    Real factor[2]={std::sqrt(.3/energy[0])-1,std::sqrt(.13/energy[1])-1};
    for(std::size_t q=0;q<fft_.size();++q){auto k=fft_.wave(q);Real k2=sq(k[0])+sq(k[1])+sq(k[2]);
        int shell=static_cast<int>(std::floor(std::sqrt(k2)+.5));FourierValue dot=0;
        for(int d=0;d<3;++d)dot+=Real(k[d])*u[d][q];
        for(int d=0;d<3;++d)u[d][q]=(shell==1 || shell==2)?factor[shell-1]*(u[d][q]-Real(k[d])*dot/k2):FourierValue{};
    }
    for(auto& a:u)fft_.inverse(a);
    std::vector<Real> slab(fft_.size()*3);
    for(std::size_t q=0;q<fft_.size();++q)for(int d=0;d<3;++d)slab[q*3+d]=u[d][q].real();
    auto delta=from_slab(slab,3);std::vector<Real> means(4);
    for(std::size_t q=0;q<cells_.size();++q){const auto& c=cells_[q];auto a=load_conservative(c.block->flow.conservative,c.index);
        means[0]+=a[0];for(int d=0;d<3;++d)means[1+d]+=a[0]*delta[q*3+d];}
    context_.mpi.sum_reals(means);Real work=0;
    for(std::size_t q=0;q<cells_.size();++q){const auto& c=cells_[q];auto b=c.index;
        auto a=load_conservative(c.block->flow.conservative,b);Real change=0;
        for(int d=0;d<3;++d){Real dv=delta[q*3+d]-(config_.hit.remove_mean_acceleration?means[d+1]/means[0]:0);
            change+=a[d+1]*dv+.5*a[0]*dv*dv;c.block->flow.conservative(b.i,b.j,b.k,d+1)+=a[0]*dv;}
        // The exact kinetic increment leaves internal energy unchanged cell by cell.
        c.block->flow.conservative(b.i,b.j,b.k,4)+=change;work+=change;
    }
    force_power_=context_.mpi.sum(work)/(means[0]*dt);
}

void HitRuntime::sample(const SimulationState& state,bool write) {
    const auto& mpi=context_.mpi;const auto& h=config_.hit;const int n=fft_.n();
    const Real count=static_cast<Real>(n)*n*n,dk=2*pi/h.length;
    const Real mu=context_.quantities.transport.viscosity(1)/context_.quantities.reference.reynolds();
    auto u=velocity_spectrum();
    std::vector<Real> shell(spectrum_integral_.size()*6),one((n/2+1)*3);
    std::vector<Real> spec(7); // energy, solenoidal, dilatational, enstrophy, div^2, helicity, E/k
    Real forced_energy=0;
    for(std::size_t q=0;q<fft_.size();++q) {
        auto k=fft_.wave(q);const Real k2=sq(k[0])+sq(k[1])+sq(k[2]),km=std::sqrt(k2);
        const auto s=static_cast<std::size_t>(std::floor(km+.5));
        Real e=0;FourierValue dot=0;
        for(int d=0;d<3;++d) {e+=.5*std::norm(u[d][q]);dot+=Real(k[d])*u[d][q];}
        const Real ed=k2>0?.5*std::norm(dot)/k2:0;
        if(k2>0) {spec[0]+=e;spec[1]+=e-ed;spec[2]+=ed;spec[6]+=e/(km*dk);}
        if(k2>0 && (h.forcing=="jhtdb_shells"?s<=2:k2<=sq(h.forcing_kmax)))forced_energy+=e-ed;
        shell[s*6]+=k2>0?e:0;shell[s*6+1]+=k2>0?e-ed:0;shell[s*6+2]+=ed;
        for(int d=0;d<3;++d) {shell[s*6+3+d]+=k2>0?.5*std::norm(u[d][q]):0;
            if(k2>0)one[static_cast<std::size_t>(std::abs(k[0]))*3+d]+=std::norm(u[d][q]);}
        // Nyquist derivative is zero for the real collocation derivative.
        for(auto& kd:k)if(std::abs(kd)==n/2)kd=0;
        std::array<FourierValue,3> omega{{FourierValue(0,dk)*(Real(k[1])*u[2][q]-Real(k[2])*u[1][q]),
          FourierValue(0,dk)*(Real(k[2])*u[0][q]-Real(k[0])*u[2][q]),
          FourierValue(0,dk)*(Real(k[0])*u[1][q]-Real(k[1])*u[0][q])}};
        dot=0;for(int d=0;d<3;++d) {spec[3]+=.5*std::norm(omega[d]);spec[5]+=(std::conj(u[d][q])*omega[d]).real();dot+=Real(k[d])*u[d][q];}
        spec[4]+=dk*dk*std::norm(dot);
    }
    mpi.sum_reals(spec);mpi.sum_reals(shell);mpi.sum_reals(one);forced_energy=mpi.sum(forced_energy);
    // Spectral derivatives are diagnostic derivatives, not the WCNS viscous stencil.
    std::vector<Real> gradients(fft_.size()*9);
    for(int d=0;d<3;++d)for(int axis=0;axis<3;++axis) {
        FourierField g=u[d];for(std::size_t q=0;q<g.size();++q) {auto k=fft_.wave(q);const int kd=std::abs(k[axis])==n/2?0:k[axis];g[q]*=FourierValue(0,dk*kd);}
        fft_.inverse(g);for(std::size_t q=0;q<g.size();++q)gradients[q*9+d*3+axis]=g[q].real();
    }
    auto grad=from_slab(gradients,9);
    // Raw volume moments; all cells have the same independently checked volume.
    std::vector<Real> raw(44);
    const int d1[6]={0,1,2,0,0,1},d2[6]={0,1,2,1,2,2};
    for(std::size_t q=0;q<cells_.size();++q) {
        const auto& c=cells_[q];const auto uc=load_conservative(c.block->flow.conservative,c.index);
        const auto primitive=temperature_primitive_from_conservative(uc,context_.quantities.gas,context_.quantities.reference,context_.quantities.floors,3);
        const Real rho=uc[0],T=primitive[4];
        const Real p=rho*T/(context_.quantities.gas.gamma()*sq(context_.quantities.reference.mach()));
        Real v[3]={uc[1]/rho,uc[2]/rho,uc[3]/rho};const Real values[6]={rho,v[0],v[1],v[2],p,T};
        for(int d=0;d<6;++d)raw[d]+=values[d];
        for(int m=0;m<6;++m){raw[6+m]+=v[d1[m]]*v[d2[m]];raw[15+m]+=rho*v[d1[m]]*v[d2[m]];}
        for(int d=0;d<3;++d)raw[12+d]+=rho*v[d];
        raw[21]+=rho*rho;raw[22]+=p*p;raw[23]+=T*T;
        const Real div=grad[q*9]+grad[q*9+4]+grad[q*9+8];Real strain=0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
            const Real sij=.5*(grad[q*9+i*3+j]+grad[q*9+j*3+i])-(i==j?div/3:0);strain+=sij*sij;
        }
        raw[24]+=2*mu*strain;raw[25]+=p*div;raw[26]+=uc[4];
        for(int d=0;d<3;++d){const Real a=grad[q*9+4*d];raw[27+3*d]+=a*a;raw[28+3*d]+=a*a*a;raw[29+3*d]+=a*a*a*a;}
    }
    mpi.sum_reals(raw);for(auto& v:raw)v/=count;
    // Centered second pass avoids catastrophic cancellation for nearly uniform
    // thermodynamic fields (and rank-count-dependent sqrt(roundoff) RMS values).
    std::vector<Real> therm_variance(3);
    for(const auto& c:cells_) {
        const auto a=temperature_primitive_from_conservative(load_conservative(c.block->flow.conservative,c.index),
            context_.quantities.gas,context_.quantities.reference,context_.quantities.floors,3);
        const Real p=a[0]*a[4]/(context_.quantities.gas.gamma()*sq(context_.quantities.reference.mach()));
        therm_variance[0]+=sq(a[0]-raw[0]);therm_variance[1]+=sq(p-raw[4]);therm_variance[2]+=sq(a[4]-raw[5]);
    }
    mpi.sum_reals(therm_variance);
    std::fill(instant_.begin(),instant_.end(),0);
    for(int d=0;d<6;++d)instant_[d]=raw[d];
    for(int m=0;m<6;++m){instant_[6+m]=raw[6+m]-raw[1+d1[m]]*raw[1+d2[m]];
        instant_[15+m]=raw[15+m]/raw[0]-raw[12+d1[m]]*raw[12+d2[m]]/sq(raw[0]);}
    for(int d=0;d<3;++d)instant_[12+d]=raw[12+d]/raw[0];
    for(int d=0;d<3;++d)instant_[21+d]=std::sqrt(therm_variance[d]/count);
    const Real K=.5*(instant_[6]+instant_[7]+instant_[8]),urms=std::sqrt(2*K/3),nu=mu/raw[0];
    const Real epsilon=raw[24]/raw[0],lambda=epsilon>0?std::sqrt(15*nu*urms*urms/epsilon):0;
    const Real eta=epsilon>0?std::pow(nu*nu*nu/epsilon,.25):0;
    instant_[24]=K;instant_[25]=.5*(raw[15]+raw[16]+raw[17])/raw[0];instant_[26]=urms;
    instant_[27]=epsilon;instant_[28]=spec[3];instant_[29]=std::sqrt(spec[4]);instant_[30]=raw[25]/raw[0];
    instant_[31]=spec[5];instant_[32]=K>0?3*pi/(4*K)*spec[6]:0;instant_[33]=lambda;
    instant_[34]=urms*lambda/nu;instant_[35]=eta;instant_[36]=n*.5*dk*eta;
    instant_[37]=std::sqrt(2*K)*context_.quantities.reference.mach()/std::sqrt(raw[5]);
    instant_[38]=spec[1];instant_[39]=spec[2];instant_[40]=forced_energy;
    instant_[41]=force_power_;instant_[42]=cooling_;instant_[43]=force_alpha_;
    for(int d=0;d<3;++d) {const Real var=raw[27+3*d];instant_[44+2*d]=var>0?raw[28+3*d]/std::pow(var,1.5):0;
        instant_[45+2*d]=var>0?raw[29+3*d]/(var*var):0;}
    instant_[50]=spec[0]-K;instant_[51]=raw[0]*std::pow(h.length,3);instant_[52]=raw[26]*std::pow(h.length,3);
    Real numerical=std::numeric_limits<Real>::quiet_NaN();
    if(sampled_ && state.time>last_sample_time_) {
        const Real dt=state.time-last_sample_time_;
        numerical=.5*(instant_[Pinindex]+previous_[Pinindex]+instant_[Pdindex]+previous_[Pdindex]
                    -instant_[Epsindex]-previous_[Epsindex])-(instant_[Kmindex]-previous_[Kmindex])/dt;
        const Real start=std::max(last_sample_time_,h.statistics_start),end=std::min(state.time,h.statistics_end);
        if(end>start && !preparing()) {
            // Linear interpolation at window edges; sample spacing is explicit in output.
            const Real a=(start-last_sample_time_)/dt,b=(end-last_sample_time_)/dt,w=end-start;
            for(std::size_t m=0;m<integral_.size();++m)integral_[m]+=w*(previous_[m]+.5*(a+b)*(instant_[m]-previous_[m]));
            // Spectrum uses right-endpoint weights; separate documented estimator.
            for(std::size_t s=0;s<spectrum_integral_.size();++s)spectrum_integral_[s]+=w*shell[s*6]/dk;
            weight_+=w;++samples_;
        }
    }
    if(write) root_io([&] {
        auto paths=output_paths();bool append=exists(paths[0]);auto f=output(paths[0],append);
        if(!append){f<<"step,time";for(const auto& name:names)f<<','<<name;f<<",epsilon_balance_residual,statistics_weight,statistics_intervals\n";}
        f<<state.step<<','<<state.time;for(Real v:instant_)f<<','<<v;f<<','<<numerical<<','<<weight_<<','<<samples_<<'\n';
        bool ap=exists(paths[1]);auto s=output(paths[1],ap);
        if(!ap)s<<"step,time,shell,k,E,E_solenoidal,E_dilatational,E_u,E_v,E_w\n";
        for(std::size_t k=0;k<shell.size()/6;++k){s<<state.step<<','<<state.time<<','<<k<<','<<k*dk;
            for(int d=0;d<6;++d)s<<','<<shell[k*6+d]/dk;
            s<<'\n';}
        bool a1=exists(paths[3]);auto l=output(paths[3],a1);
        if(!a1)l<<"step,time,kx,F11,F22,F33\n";
        for(int k=0;k<=n/2;++k){l<<state.step<<','<<state.time<<','<<k*dk;for(int d=0;d<3;++d)l<<','<<one[k*3+d]/dk;l<<'\n';}
        bool ac=exists(paths[4]);auto cr=output(paths[4],ac);
        if(!ac)cr<<"step,time,r,R11,R22,R33,D11,D22,D33\n";
        for(int r=0;r<=n/2;++r){Real corr[3]={};for(int k=0;k<=n/2;++k)for(int d=0;d<3;++d)corr[d]+=one[k*3+d]*std::cos(2*pi*k*r/n);
            cr<<state.step<<','<<state.time<<','<<r*h.length/n;for(Real c:corr)cr<<','<<c;
            for(int d=0;d<3;++d)cr<<','<<2*(instant_[6+d]-corr[d]);
            cr<<'\n';}
    });
    previous_=instant_;sampled_=true;last_sample_time_=state.time;last_sample_step_=state.step;
}
void HitRuntime::root_io(const std::function<void()>& action) const {
    std::string error;if(context_.mpi.rank()==0)try{action();}catch(const std::exception& e){error=e.what();}
    error=context_.mpi.broadcast_string(error);if(!error.empty())throw std::runtime_error(error);
}
std::vector<std::string> HitRuntime::output_paths() const {
    std::vector<std::string> p;for(const auto* name:{"history","spectrum","means","spectrum_1d","correlation","spectrum_mean","metadata"})
        p.push_back(config_.output.directory+(preparing()?"/preparation/":"/")+config_.case_name+"_hit_"+name+(std::string(name)=="metadata"?".txt":".csv"));
    return p;
}
void HitRuntime::write_means() const {
    root_io([&]{auto paths=output_paths();auto f=output(paths[2]);f<<"quantity,time_mean,weight,intervals\n";
        for(std::size_t m=0;m<integral_.size();++m)f<<names[m]<<','<<(weight_>0?integral_[m]/weight_:std::numeric_limits<Real>::quiet_NaN())<<','<<weight_<<','<<samples_<<'\n';
        auto s=output(paths[5]);s<<"shell,k,E_time_mean,weight\n";
        for(std::size_t k=0;k<spectrum_integral_.size();++k)s<<k<<','<<k*2*pi/config_.hit.length<<','
            <<(weight_>0?spectrum_integral_[k]/weight_:std::numeric_limits<Real>::quiet_NaN())<<','<<weight_<<'\n';});
}
void HitRuntime::on_initial(const SimulationState& s) {
    if(restored_ && (last_step_!=s.step || std::abs(last_step_time_-s.time)>1e-12))throw std::runtime_error("HIT checkpoint time mismatch");
    last_step_=s.step;last_step_time_=s.time;
    if(!restored_)sample(s,true);
    root_io([&]{auto f=output(output_paths()[6]);f<<"WCNS v2.6\n"<<config_.hit.signature()<<"\nspectrum_content_fnv1a="<<spectrum_identity_
        <<"\nFFT="<<fft_.backend()<<"\nMPI_ranks="<<context_.mpi.size()
        <<"\ninitialization="<<config_.hit.initialization<<"\nphase="<<(preparing()?"preparation":"production")
        <<"\npreparation_complete="<<preparation_complete_<<"\npreparation_time="<<config_.hit.preparation_time
        <<"\npreparation_steps="<<preparation_steps_
        <<"\nforward=DFT/N^3; inverse=sum; shell=round(|integer wave|); E=shell kinetic energy/dk\n"
        <<"gradients=spectral real collocation; Nyquist derivative=0\n"
        <<"statistics=volume ensemble then trapezoidal time integration; mean spectrum=right endpoint\n"
        <<"epsilon_balance_residual includes spatial and temporal discretization and sampled budget error\n";
        if(config_.hit.initialization=="analytic_random_phase") {
            const auto& h=config_.hit;
            f<<"analytic_A="<<h.spectrum_amplitude<<"\nanalytic_k0="<<h.peak_wave
             <<"\nanalytic_K0="<<3*h.spectrum_amplitude/64*std::sqrt(2*pi)*std::pow(h.peak_wave,5)
             <<"\nanalytic_tau="<<std::sqrt(32/h.spectrum_amplitude)*std::pow(2*pi,.25)*std::pow(h.peak_wave,-3.5)<<'\n';
        }
    });
}
void HitRuntime::on_step(const SimulationState& s,bool) {
    if(!s.residuals.finite || s.step<=last_step_ || s.time<=last_step_time_)return;
    bool due=s.step%config_.hit.sample_every_steps==0;
    for(Real t:config_.hit.sample_times)due=due || (last_step_time_<t && s.time>=t-1e-12*std::max(Real{1},t));
    last_step_=s.step;last_step_time_=s.time;
    if(due){sample(s,true);if(samples_%config_.hit.write_every_samples==0)write_means();}
}
void HitRuntime::on_final(const SimulationState& s) {
    if(s.stop_reason==StopReason::NumericalFailure)return;
    if(!sampled_ || s.step!=last_sample_step_)sample(s,true);
    write_means();
}
Real HitRuntime::next_time_event(const SimulationState& s) const {
    if(preparing())return config_.hit.preparation_time;
    for(Real t:config_.hit.sample_times)if(t>s.time+1e-12*std::max(Real{1},s.time))return t;
    return std::numeric_limits<Real>::infinity();
}
std::string HitRuntime::serialize() const {
    std::ostringstream s;s<<std::setprecision(17)<<config_.hit.signature()<<'\n'<<spectrum_identity_<<'\n'
      <<weight_<<' '<<samples_<<' '<<last_sample_time_<<' '<<last_sample_step_<<' '<<last_step_time_<<' '<<last_step_<<' '<<sampled_<<' '<<force_power_<<'\n';
    for(const auto* a:{&integral_,&previous_,&spectrum_integral_}){s<<a->size()<<' ';for(Real v:*a)s<<v<<' ';s<<'\n';}
    if(config_.hit.preparation_time>0)s<<"preparation_v1 "<<preparation_complete_<<' '<<preparation_steps_<<'\n';
    return s.str();
}
void HitRuntime::restore(const std::string& text) {
    std::istringstream s(text);std::string signature,identity;std::getline(s,signature);std::getline(s,identity);
    Real w=0,ts=0,t=0;std::size_t samples=0,ss=0,step=0;bool sampled=false;Real power=0;s>>w>>samples>>ts>>ss>>t>>step>>sampled>>power;
    if(!s || !std::isfinite(power) || signature!=config_.hit.signature() || identity!=spectrum_identity_ || !std::isfinite(w) || w<0
       || !std::isfinite(ts) || ts<0 || !std::isfinite(t) || t<ts || step<ss)throw std::runtime_error("invalid HIT checkpoint header or spectrum content");
    auto a=integral_,b=previous_,c=spectrum_integral_;
    for(auto* v:{&a,&b,&c}){std::size_t n=0;s>>n;if(n!=v->size())throw std::runtime_error("HIT checkpoint dimensions");
        for(auto& x:*v){s>>x;if(!s || !std::isfinite(x))throw std::runtime_error("invalid HIT checkpoint statistics");}}
    bool completed=true;std::size_t prepared_steps=0;
    if(config_.hit.preparation_time>0) {
        std::string tag;int flag=-1;s>>tag>>flag>>prepared_steps;
        if(!s || tag!="preparation_v1" || (flag!=0 && flag!=1)
           || (!flag && (prepared_steps!=0 || t>config_.hit.preparation_time+1e-12))
           || (flag && prepared_steps==0))throw std::runtime_error("invalid HIT preparation checkpoint state");
        completed=flag==1;
    }
    std::string extra;if(s>>extra)throw std::runtime_error("trailing HIT checkpoint state");
    integral_=std::move(a);previous_=std::move(b);spectrum_integral_=std::move(c);
    weight_=w;samples_=samples;last_sample_time_=ts;last_sample_step_=ss;last_step_time_=t;last_step_=step;sampled_=sampled;force_power_=power;restored_=true;
    preparation_complete_=completed;preparation_steps_=prepared_steps;
}
} // namespace wcns
