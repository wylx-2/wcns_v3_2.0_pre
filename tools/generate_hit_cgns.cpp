#include <cgnslib.h>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
void check(int status){if(status!=CG_OK)throw std::runtime_error(cg_get_error());}
int integer(const char* text){std::size_t used=0;int n=std::stoi(text,&used);
 if(used!=std::string(text).size() || n<1 || n>512)throw std::invalid_argument("invalid grid integer");
 return n;}
std::string name(int x,int y,int z,int b){return "HIT_"+std::to_string((z*b+y)*b+x);}
}
int main(int argc,char** argv){int file=0;
 try {
  if(argc!=5)throw std::invalid_argument("usage: wcns_generate_hit_cgns output.cgns N blocks_per_axis length");
  if(std::ifstream(argv[1]).good())throw std::invalid_argument("refusing to overwrite mesh");
  int n=integer(argv[2]),b=integer(argv[3]);std::size_t used=0;double L=std::stod(argv[4],&used);
  if(n<8 || (n&(n-1)) || n%b || n/b<4 || !std::isfinite(L) || L<=0 || used!=std::string(argv[4]).size())
    throw std::invalid_argument("require power-of-two N, equal blocks with >=4 cells and positive length");
  const int c=n/b,v=c+1;check(cg_set_file_type(CG_FILE_ADF));check(cg_open(argv[1],CG_MODE_WRITE,&file));
  int base=0;check(cg_base_write(file,"HIT",3,3,&base));
  for(int z=0;z<b;++z)for(int y=0;y<b;++y)for(int x=0;x<b;++x){
   int zone=0,coord=0;cgsize_t sizes[9]={v,v,v,c,c,c,0,0,0};
   check(cg_zone_write(file,base,name(x,y,z,b).c_str(),sizes,Structured,&zone));
   std::vector<double> data(static_cast<std::size_t>(v)*v*v);int origin[3]={x*c,y*c,z*c};
   const char* coordinates[3]={"CoordinateX","CoordinateY","CoordinateZ"};
   for(int d=0;d<3;++d){for(int k=0;k<v;++k)for(int j=0;j<v;++j)for(int i=0;i<v;++i){int a[3]={i,j,k};
      data[(static_cast<std::size_t>(k)*v+j)*v+i]=L*(origin[d]+a[d])/n;}
     check(cg_coord_write(file,base,zone,RealDouble,coordinates[d],data.data(),&coord));}
   for(int d=0;d<3;++d)for(int side=0;side<2;++side){
    int donor[3]={x,y,z};donor[d]+=side?1:-1;bool periodic=donor[d]<0 || donor[d]>=b;donor[d]=(donor[d]+b)%b;
    std::array<cgsize_t,6> r{{1,1,1,v,v,v}},dr=r;r[d]=r[d+3]=side?v:1;dr[d]=dr[d+3]=side?1:v;
    int transform[3]={1,2,3},id=0;std::string connection="axis"+std::to_string(d)+(side?"plus":"minus");
    check(cg_1to1_write(file,base,zone,connection.c_str(),name(donor[0],donor[1],donor[2],b).c_str(),r.data(),dr.data(),transform,&id));
    if(periodic){float zero[3]={},shift[3]={};shift[d]=static_cast<float>(side?-L:L);
      check(cg_1to1_periodic_write(file,base,zone,id,zero,zero,shift));}
   }
  }
  check(cg_close(file));file=0;
  std::cout<<std::setprecision(17)<<"HIT mesh N="<<n<<" blocks="<<b*b*b<<" length="<<L<<" cells="<<static_cast<std::size_t>(n)*n*n<<'\n';return 0;
 }catch(const std::exception& e){if(file)cg_close(file);std::cerr<<e.what()<<'\n';return 1;}
}
