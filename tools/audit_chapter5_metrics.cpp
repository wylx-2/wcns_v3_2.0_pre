#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <iostream>
#include <iomanip>
int main(int argc,char**argv) {try {
 if(argc<2 || argc>3) return 2;
 wcns::CgnsReader reader;const auto metadata=reader.read_metadata(argv[1]);
 wcns::MetricBuildOptions options;
 if(argc==3) options.maximum_reference_relative_difference=std::stod(argv[2]);
 const auto profile=wcns::ProfileFactory::from_string("scmm6_wcns");
 std::size_t cells=0;double maxdiff=0,minjac=1e300;
 for(const auto& z:metadata.zones) {
    auto block=reader.read_block(argv[1],z,0,0);
    auto m=wcns::initialize_metric_field(block,profile,options);
    const auto e=block.cell_extent();cells+=static_cast<std::size_t>(e.ni)*e.nj*e.nk;
    for(int k=0;k<e.nk;++k) for(int j=0;j<e.nj;++j) for(int i=0;i<e.ni;++i) minjac=std::min(minjac,m.metric.jacobian()(i,j,k));
    maxdiff=std::max(maxdiff,m.diagnostics.maximum_jacobian_relative_difference);
    std::cout<<std::setprecision(17)<<"zone="<<z.name<<" maximum_relative_difference="<<m.diagnostics.maximum_jacobian_relative_difference<<" fallback="<<m.diagnostics.fallback_cell_count<<std::endl;
 }
 std::cout<<"PASS: zones="<<metadata.zones.size()<<" cells="<<cells<<" min_scmm_jacobian="<<minjac<<" max_reference_difference="<<maxdiff<<'\n';return 0;
 }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
