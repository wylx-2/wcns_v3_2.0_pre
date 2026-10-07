#include "test_support.hpp"
#include <wcns/physics/periodic_hill.hpp>
#include <wcns/runtime/periodic_hill.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

int main(int argc,char** argv) {
    using namespace wcns;
    try {
        MpiRuntime mpi(argc,argv);
        WCNS_REQUIRE(mpi.size()==1);
        WCNS_REQUIRE_NEAR(periodic_hill_geometry(0)[0],1,1e-14);
        WCNS_REQUIRE_NEAR(periodic_hill_geometry(4.5)[0],0,1e-14);
        for(Real x:{.4,.8,1.2,1.8,3.0}) {
            WCNS_REQUIRE_NEAR(periodic_hill_geometry(x)[0],periodic_hill_geometry(9-x)[0],1e-12);
            const Real d=(periodic_hill_geometry(x+1e-6)[0]-periodic_hill_geometry(x-1e-6)[0])/2e-6;
            WCNS_REQUIRE_NEAR(d,periodic_hill_geometry(x)[1],1e-7);
        }
        WCNS_REQUIRE_NEAR(hill_wall_derivative(.3+.02,.6+.08,.1,.2),3,1e-13);
        WCNS_REQUIRE_NEAR(hill_feedback_force(.01,.2,.1,.05,1,1,.5),.12,1e-14);
        WCNS_REQUIRE_NEAR(hill_feedback_force(.49,2,0,.1,1,1,.5),.5,0);
        CaseConfig cfg; cfg.periodic_hill.enabled=true; cfg.periodic_hill.statistics_start=0;
        cfg.periodic_hill.statistics_end=10; cfg.periodic_hill.write_every_steps=100;
        cfg.output.directory="."; cfg.case_name="unit-hill";
        GasModelInput gi; gi.specific_gas_constant=1;
        const auto gas=GasModel::from_input(gi);
        ReferenceInput ri; ri.velocity=1; ri.density=1; ri.length=1; ri.temperature=1/1.4/.01; ri.viscosity=1./1400;
        const auto ref=ReferenceScales::derive(ri,gas);
        InitialConditionConfig initial; initial.type="periodic_hill";
        initial.parameters={{"bulk_velocity",1},{"perturbation_amplitude",.05}};
        for(Real x:{0.,.4,1.5,4.5,8.8,9.}) for(Real y:{periodic_hill_geometry(x)[0],3.035}) {
            const auto p=FlowInitializer::evaluate(initial,{x,y,1},gas,ref,{},3);
            for(int a=1;a<=3;++a) WCNS_REQUIRE_NEAR(p[a],0,1e-14);
        }
        SourceTermConfig src; src.enable_source_terms=true; src.models={SourceModelKind::PressureGradient};
        auto registry=SourceTermRegistry::create_stage_j(src); registry.set_pressure_gradient_x(.25);
        const auto force=registry.evaluate({2,6,0,0,100},{0,0,0},0,3);
        WCNS_REQUIRE_NEAR(force[1],.25,0); WCNS_REQUIRE_NEAR(force[4],.75,0);
        const int nx=32,ny=12,nz=8;
        StructuredBlock b(0,"hill",0,3,3,{nx+1,ny+1,nz+1},3);
        for(int k=0;k<=nz;++k) for(int j=0;j<=ny;++j) for(int i=0;i<=nx;++i) {
            const Real x=9.*i/nx,wall=periodic_hill_geometry(x)[0];
            b.coordinates.x(i,j,k)=x; b.coordinates.y(i,j,k)=wall+(3.035-wall)*j/ny; b.coordinates.z(i,j,k)=4.5*k/nz;
        }
        const auto profile=ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
        auto metric=initialize_metric_field(b,profile).metric;
        PartitionConfig pc; const auto plan=StructuredPartitionPlan::build({{0,"hill",3,{nx,ny,nz}}},1,pc);
        LocalBlockSet local(0,{std::move(b)},plan.distribution());
        BlockMetricMap metrics; metrics.emplace(0,std::move(metric));
        GlobalConservationWeights unused;
        StatisticContext ctx{mpi,local,metrics,plan,unused,profile,{gas,ref,{},TransportModel({}),false},nullptr,true};
        auto fill=[&](Real base,Real amplitude) {
            auto& block=local.block(0);
            for(int k=0;k<nz;++k) for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) {
                // Include a spatial mean gradient: it must not create profile variance.
                const Real x=metrics.at(0).cell_coordinates().x(i,j,k);
                const Real u=base+.01*x+amplitude*(k%2 ? 1:-1);
                store_state(block.flow.conservative,{i,j,k},thermodynamic_conservative({1,u,0,0,1},gas,ref,{},3));
            }
        };
        fill(2,1); PeriodicHillRuntime hill(cfg,ctx); SimulationState s; hill.on_initial(s);
        s.step=1; s.time=.1; s.residuals.finite=true; hill.on_step(s,false);
        const auto checkpoint=hill.serialize(); PeriodicHillRuntime resumed(cfg,ctx); resumed.restore(checkpoint); resumed.on_initial(s);
        fill(4,2); s.step=2; s.time=.4; hill.on_step(s,false); resumed.on_step(s,false);
        WCNS_REQUIRE(hill.serialize()==resumed.serialize());
        WCNS_REQUIRE_NEAR(hill.weight(),.4,1e-14);
        resumed.on_step(s,false); WCNS_REQUIRE_NEAR(resumed.weight(),.4,1e-14);
        WCNS_REQUIRE_THROWS(std::runtime_error,resumed.restore(checkpoint+"garbage"));
        hill.on_final(s);
        for(const auto& path:{hill.output_paths()[1],hill.output_paths()[2]}) {
            std::ifstream f(path); std::string line; std::getline(f,line); int rows=0;
            while(std::getline(f,line)) {
                std::replace(line.begin(),line.end(),',',' '); std::istringstream r(line); Real v[14];
                for(auto& x:v) r>>x; WCNS_REQUIRE(bool(r));
                if(path==hill.output_paths()[1]) WCNS_REQUIRE_NEAR(v[3],3.5+.01*v[0],1e-12);
                WCNS_REQUIRE_NEAR(v[8],4,1e-11); ++rows;
            }
            WCNS_REQUIRE(rows>0);
        }
        cfg.periodic_hill.statistics_start=.05; cfg.periodic_hill.statistics_end=.3;
        PeriodicHillRuntime clipped(cfg,ctx); s.step=0;s.time=0;clipped.on_initial(s);
        s.step=1;s.time=.1;clipped.on_step(s,false);s.step=2;s.time=.4;clipped.on_step(s,false);
        WCNS_REQUIRE_NEAR(clipped.weight(),.25,1e-14);
        std::cout<<"periodic hill geometry, initialization, source work, PI, span/time moments, window and restart: PASS\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
