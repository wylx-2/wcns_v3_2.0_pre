#include "test_support.hpp"
#include <wcns/runtime/hit.hpp>
#include <wcns/solver/euler.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <cstdio>

int main(int argc,char** argv){try {
    using namespace wcns;MpiRuntime mpi(argc,argv);WCNS_REQUIRE(mpi.size()==1);
    const int n=8;const Real L=2*std::acos(-1.);
    GasModelInput gi;gi.specific_gas_constant=1;auto gas=GasModel::from_input(gi);
    ReferenceInput ri;ri.velocity=1;ri.density=1;ri.temperature=1/(1.4*.01);ri.length=1;ri.viscosity=.001;
    auto ref=ReferenceScales::derive(ri,gas);auto profile=ProfileFactory::from_string("scmm6_wcns");
    auto plan=StructuredPartitionPlan::build({{0,"cube",3,{n,n,n}}},1,PartitionConfig{});
    StructuredBlock b(0,"cube",0,3,3,{n+1,n+1,n+1},3);
    for(int k=0;k<=n;++k)for(int j=0;j<=n;++j)for(int i=0;i<=n;++i){b.coordinates.x(i,j,k)=L*i/n;b.coordinates.y(i,j,k)=L*j/n;b.coordinates.z(i,j,k)=L*k/n;}
    BlockMetricMap metrics;metrics.emplace(0,initialize_metric_field(b,profile).metric);
    std::vector<StructuredBlock> blocks;blocks.push_back(std::move(b));LocalBlockSet local(0,std::move(blocks),plan.distribution());
    GlobalConservationWeights unused;StatisticContext ctx{mpi,local,metrics,plan,unused,profile,{gas,ref,{},TransportModel({}),false},nullptr,true};
    CaseConfig cfg;cfg.hit.type="forced";cfg.hit.n=n;cfg.hit.length=L;cfg.hit.sample_every_steps=1;cfg.hit.statistics_end=2;
    cfg.case_name="unit-hit";cfg.output.directory=".";
    auto& block=local.block(0);const auto& xyz=metrics.at(0).cell_coordinates();
    auto fill=[&]{for(int k=0;k<n;++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
        Real x=xyz.x(i,j,k),y=xyz.y(i,j,k),z=xyz.z(i,j,k);
        store_state(block.flow.conservative,{i,j,k},thermodynamic_conservative({1,std::sin(y)+.5*std::sin(2*y),std::cos(z)+.5*std::cos(2*z),std::sin(x)+.5*std::sin(2*x),1},gas,ref,{},3));
        for(int d=0;d<5;++d)block.flow.residual(i,j,k,d)=0;
    }};
    fill();HitRuntime runtime(cfg,ctx,local,"");for(auto& p:runtime.output_paths())std::remove(p.c_str());runtime.add_stage_source(0);
    Real power=0,momentum=0;
    for(int k=0;k<n;++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){power+=block.flow.residual(i,j,k,4);momentum+=block.flow.residual(i,j,k,1);}
    WCNS_REQUIRE_NEAR(power/(n*n*n),.103,1e-13);WCNS_REQUIRE_NEAR(momentum,0,1e-12);
    SimulationState s;runtime.on_initial(s);s.step=1;s.time=.5;s.residuals.finite=true;runtime.on_step(s,false);runtime.on_final(s);
    std::ifstream f(runtime.output_paths()[0]);std::string header,line;std::getline(f,header);std::getline(f,line);
    std::istringstream h(header),v(line);std::string key,value;std::map<std::string,Real> row;
    while(std::getline(h,key,',') && std::getline(v,value,','))row[key]=std::stod(value);
    WCNS_REQUIRE_NEAR(row.at("K"),.9375,1e-13);WCNS_REQUIRE_NEAR(row.at("epsilon_viscous"),.003,1e-14);
    WCNS_REQUIRE_NEAR(row.at("enstrophy"),1.5,1e-13);WCNS_REQUIRE_NEAR(row.at("divergence_rms"),0,1e-12);
    WCNS_REQUIRE_NEAR(row.at("uu"),.625,1e-13);WCNS_REQUIRE_NEAR(row.at("p_rms"),0,1e-10);
    HitRuntime restored(cfg,ctx,local,"");restored.restore(runtime.serialize());restored.on_initial(s);
    WCNS_REQUIRE(runtime.serialize()==restored.serialize());WCNS_REQUIRE_THROWS(std::runtime_error,restored.restore(runtime.serialize()+"bad"));
    f.close();for(auto& p:runtime.output_paths())std::remove(p.c_str());
    cfg.hit.forcing="jhtdb_shells";HitRuntime shells(cfg,ctx,local,"");fill();shells.accepted_step_transform(.001);
    Real K=0;for(int k=0;k<n;++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
        const auto u=load_conservative(block.flow.conservative,{i,j,k});for(int d=1;d<4;++d)K+=.5*u[d]*u[d]/u[0];
        auto p=temperature_primitive_from_conservative(u,gas,ref,{},3);WCNS_REQUIRE_NEAR(p[4],1,1e-13);}
    WCNS_REQUIRE_NEAR(K/(n*n*n),.43,2e-13);
    cfg.hit.type="decay";cfg.hit.initialization="analytic_random_phase";cfg.hit.peak_wave=1.5;
    cfg.hit.validate();HitRuntime analytic(cfg,ctx,local,"");analytic.initialize();
    Real initial_K=0;std::array<Real,3> mean{};
    for(int k=0;k<n;++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
        const auto u=load_conservative(block.flow.conservative,{i,j,k});
        for(int d=0;d<3;++d){initial_K+=.5*u[d+1]*u[d+1];mean[d]+=u[d+1];}
        const auto p=temperature_primitive_from_conservative(u,gas,ref,{},3);
        WCNS_REQUIRE_NEAR(p[4],1,1e-13);WCNS_REQUIRE_NEAR(u[0],1,1e-15);
    }
    const Real target_K=3*cfg.hit.spectrum_amplitude/64*std::sqrt(L)*std::pow(cfg.hit.peak_wave,5);
    WCNS_REQUIRE_NEAR(initial_K/(n*n*n),target_K,1e-14);
    for(Real x:mean)WCNS_REQUIRE_NEAR(x/(n*n*n),0,1e-14);
    HitConfig invalid_preparation=cfg.hit;invalid_preparation.preparation_time=.01;
    WCNS_REQUIRE_THROWS(CaseConfigurationError,invalid_preparation.validate());
    HitConfig bad=cfg.hit;bad.n=17;WCNS_REQUIRE_THROWS(CaseConfigurationError,bad.validate());
    std::cout<<"HIT analytic moments, gradients, source work, exact shell rescaling and checkpoint passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
