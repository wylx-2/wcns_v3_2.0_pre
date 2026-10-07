#include "test_support.hpp"
#include <wcns/runtime/chapter5.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/physics/chapter5.hpp>
#include <wcns/solver/viscous_boundary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>

int main(int argc,char**argv) {try {
 using namespace wcns;MpiRuntime mpi(argc,argv);WCNS_REQUIRE(mpi.size()==1);
 for(int n=1;n<argc;++n) {
    const auto config=CaseConfig::from_file(argv[n]);config.validate();
    WCNS_REQUIRE(config.chapter5.enabled());
    std::cout<<"validated configuration "<<argv[n]<<'\n';
 }
 WCNS_REQUIRE_NEAR(ramp_inlet(0)[1],0,1e-12);WCNS_REQUIRE_NEAR(ramp_inlet(0)[4],1.9,1e-12);
 WCNS_REQUIRE_NEAR(ramp_inlet(1)[1],.99,2e-7);WCNS_REQUIRE_NEAR(ramp_inlet(4)[1],1,1e-12);
 for(Real y=0;y<=4;y+=.01) {auto a=ramp_inlet(y);WCNS_REQUIRE_NEAR(a[0]*a[4],1,2e-15);WCNS_REQUIRE(a[1]>=-1e-12 && a[1]<=1+1e-12);}
 SourceTermConfig st;st.enable_source_terms=true;st.models={SourceModelKind::RampTrip};st.ramp_trip_amplitude=5;
 auto source=SourceTermRegistry::create_stage_j(st);
 auto a=source.evaluate({2,6,0,0,100},{2.5,.005,.25},.5,3);
 WCNS_REQUIRE(a[1]<0);WCNS_REQUIRE_NEAR(a[4],3*a[1],1e-14);WCNS_REQUIRE_NEAR(a[0],0,0);
 WCNS_REQUIRE_NEAR(ramp_trip(2.5,.005,.25,.5,5,6),ramp_trip(2.5,.005,6.25,.5,5,6),1e-13);
 WCNS_REQUIRE_NEAR(ramp_trip(2.2,.005,.25,.5,5,6),0,0);
 WCNS_REQUIRE_NEAR(ramp_trip(2.5,.02,.25,.5,5,6),0,0);
 GasModelInput gi;gi.specific_gas_constant=1;auto gas=GasModel::from_input(gi);
 ReferenceInput ri;ri.velocity=1;ri.density=1;ri.temperature=1/1.4/.04;ri.length=1;ri.viscosity=1./60000;
 auto ref=ReferenceScales::derive(ri,gas);const auto profile=ProfileFactory::from_string("scmm6_wcns");
 // Regression: SCMM cell centres above a concave parabola used to be paired
 // with a straight-chord centroid, making their apparent wall distance negative.
 StructuredBlock curved(0,"curved",0,2,2,{9,9,1},3);
 for(int j=0;j<=8;++j) for(int i=0;i<=8;++i) {
    Real x=.01*i;curved.coordinates.x(i,j,0)=x;curved.coordinates.y(i,j,0)=x*x+.000001*j;curved.coordinates.z(i,j,0)=0;
 }
 BoundaryPatch patch;patch.name="wall";patch.face={Axis::J,Side::Lower};patch.type=BoundaryType::NoSlipAdiabaticWall;
 auto cm=initialize_metric_field(curved,profile).metric;
 const auto point=viscous_wall_coordinates(curved,cm,patch,{3,0,0},profile);
 WCNS_REQUIRE_NEAR(point[0],.035,1e-12);WCNS_REQUIRE_NEAR(point[1],.035*.035,1e-12);
 WCNS_REQUIRE_NEAR(cm.cell_coordinates().y(3,0,0)-point[1],.5e-6,1e-12);
 std::vector<PartitionZone> pz{{0,"one",3,{8,8,8}},{1,"two",3,{8,8,8}}};
 PartitionConfig pc;auto plan=StructuredPartitionPlan::build(pz,1,pc);
 std::vector<StructuredBlock> blocks;BlockMetricMap metrics;
 for(int z=0;z<2;++z) {
    StructuredBlock b(z,z?"two":"one",0,3,3,{9,9,9},3);
    for(int k=0;k<=8;++k) for(int j=0;j<=8;++j) for(int i=0;i<=8;++i) {
        b.coordinates.x(i,j,k)=z+i/8.;b.coordinates.y(i,j,k)=j/8.;b.coordinates.z(i,j,k)=.2*k/8.;
    }
    b.boundaries.push_back(patch);metrics.emplace(z,initialize_metric_field(b,profile).metric);blocks.push_back(std::move(b));
 }
 LocalBlockSet local(0,std::move(blocks),plan.distribution());GlobalConservationWeights unused;
 StatisticContext ctx{mpi,local,metrics,plan,unused,profile,{gas,ref,{},TransportModel({}),false},nullptr,true};
 CaseConfig cfg;cfg.chapter5.type="sd7003";cfg.chapter5.statistics_start=0;cfg.chapter5.statistics_end=4;
 cfg.chapter5.write_every_steps=100;cfg.chapter5.history_every_steps=100;cfg.case_name="unit-chapter5";cfg.output.directory=".";
 auto fill=[&](Real shift) {for(int z=0;z<2;++z) {auto& b=local.block(z);for(int k=0;k<8;++k) for(int j=0;j<8;++j) for(int i=0;i<8;++i) {
    const Real r=k%2?3:1,u=r+shift;store_state(b.flow.conservative,{i,j,k},thermodynamic_conservative({r,u,0,0,1},gas,ref,{},3));
 }}};
 fill(0);Chapter5Runtime run(cfg,ctx);SimulationState state;run.on_initial(state);
 state.step=1;state.time=1;state.residuals.finite=true;run.on_step(state,false);
 Chapter5Runtime restart(cfg,ctx);restart.restore(run.serialize());restart.on_initial(state);
 fill(2);state.step=2;state.time=4;run.on_step(state,false);restart.on_step(state,false);
 WCNS_REQUIRE(run.serialize()==restart.serialize());const auto serialized=run.serialize();run.on_step(state,false);WCNS_REQUIRE(run.serialize()==serialized);
 WCNS_REQUIRE_THROWS(std::runtime_error,restart.restore(serialized+"garbage"));
 run.on_final(state);
 std::ifstream csv(run.output_paths()[2]);std::string line;std::getline(csv,line);std::getline(csv,line);std::istringstream row(line);std::vector<double> values;
 while(std::getline(row,line,',')) values.push_back(std::stod(line));
 WCNS_REQUIRE_NEAR(values[11],2,1e-13);WCNS_REQUIRE_NEAR(values[12],3.5,1e-13);
 WCNS_REQUIRE_NEAR(values[17],1.75,1e-12);WCNS_REQUIRE_NEAR(values[23],4,1e-13);WCNS_REQUIRE_NEAR(values[26],1.5,1e-12);
 WCNS_REQUIRE_NEAR(values[32],1,1e-13);WCNS_REQUIRE_NEAR(values[35],4,1e-13);
 for(const auto& path:run.output_paths()) std::remove(path.c_str());
 // Exact gradient on a skew grid: u=-2y,v=2x gives omega_z=4,Q=4.
 auto& b=local.block(0);const auto& c=metrics.at(0).cell_coordinates();
 for(int k=0;k<8;++k) for(int j=0;j<8;++j) for(int i=0;i<8;++i)
    store_state(b.flow.conservative,{i,j,k},thermodynamic_conservative({1,-2*c.y(i,j,k),2*c.x(i,j,k),0,1},gas,ref,{},3));
 auto registry=FieldQuantityRegistry::create_builtin();
 auto q=registry.evaluate("Q",b,metrics.at(0),ctx.quantities);auto omega=registry.evaluate("vorticity_z",b,metrics.at(0),ctx.quantities);
 // QuantityField stores the same cell-shaped values used by the field writer.
 WCNS_REQUIRE_NEAR(q.values[(3*8+3)*8+3],4,1e-12);WCNS_REQUIRE_NEAR(omega.values[(3*8+3)*8+3],4,1e-12);
 std::cout<<"chapter5 inlet/source/curved-wall/multizone moments/restart/Q passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
