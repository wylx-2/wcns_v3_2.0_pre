#include <wcns/physics/periodic_hill.hpp>
#include <cgnslib.h>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(int status) { if(status!=CG_OK) throw std::runtime_error(cg_get_error()); }
int count(const char* s) {
    std::size_t n=0; const int v=std::stoi(s,&n);
    if(n!=std::string(s).size() || v<8 || v>2048) throw std::invalid_argument("cell counts must be in [8,2048]");
    return v;
}
}
int main(int argc,char** argv) {
    int file=0;
    try {
        if(argc!=2 && argc!=5) throw std::invalid_argument("usage: wcns_generate_periodic_hill_cgns output.cgns [nx ny nz]; defaults 88 48 24");
        if(std::ifstream(argv[1]).good()) throw std::invalid_argument("refusing to overwrite existing mesh");
        const int nx=argc==5 ? count(argv[2]) : 88, ny=argc==5 ? count(argv[3]) : 48, nz=argc==5 ? count(argv[4]) : 24;
        const std::size_t total=static_cast<std::size_t>(nx+1)*(ny+1)*(nz+1);
        if(total>50000000) throw std::invalid_argument("mesh exceeds generator memory limit");
        const double pi=std::acos(-1.0),beta=2.2,cluster=.45;
        std::vector<double> x(total),y(total),z(total);
        for(int k=0;k<=nz;++k) for(int j=0;j<=ny;++j) for(int i=0;i<=nx;++i) {
            const double xi=double(i)/nx,eta=double(j)/ny;
            const double X=9*(xi-cluster*std::sin(2*pi*xi)/(2*pi));
            const double f=.5*(1+std::tanh(beta*(2*eta-1))/std::tanh(beta));
            const double wall=wcns::periodic_hill_geometry(X)[0];
            const auto q=(static_cast<std::size_t>(k)*(ny+1)+j)*(nx+1)+i;
            x[q]=X; y[q]=wall+(3.035-wall)*f; z[q]=4.5*k/nz;
        }
        check(cg_set_file_type(CG_FILE_ADF)); check(cg_open(argv[1],CG_MODE_WRITE,&file));
        int base=0,zone=0,coord=0; check(cg_base_write(file,"Base",3,3,&base));
        cgsize_t size[9]={nx+1,ny+1,nz+1,nx,ny,nz,0,0,0};
        check(cg_zone_write(file,base,"PeriodicHill",size,Structured,&zone));
        check(cg_coord_write(file,base,zone,RealDouble,"CoordinateX",x.data(),&coord));
        check(cg_coord_write(file,base,zone,RealDouble,"CoordinateY",y.data(),&coord));
        check(cg_coord_write(file,base,zone,RealDouble,"CoordinateZ",z.data(),&coord));
        for(int side=0;side<2;++side) {
            cgsize_t range[6]={1,side ? ny+1:1,1,nx+1,side ? ny+1:1,nz+1}; int bc=0;
            check(cg_boco_write(file,base,zone,side ? "top":"bottom",BCWallViscous,PointRange,2,range,&bc));
        }
        auto connect=[&](const char* name,std::array<cgsize_t,6> range,std::array<cgsize_t,6> donor,std::array<float,3> shift) {
            int transform[3]={1,2,3},id=0; float zero[3]={0,0,0};
            check(cg_1to1_write(file,base,zone,name,"PeriodicHill",range.data(),donor.data(),transform,&id));
            check(cg_1to1_periodic_write(file,base,zone,id,zero,zero,shift.data()));
        };
        connect("imin",{1,1,1,1,ny+1,nz+1},{nx+1,1,1,nx+1,ny+1,nz+1},{9,0,0});
        connect("imax",{nx+1,1,1,nx+1,ny+1,nz+1},{1,1,1,1,ny+1,nz+1},{-9,0,0});
        connect("kmin",{1,1,1,nx+1,ny+1,1},{1,1,nz+1,nx+1,ny+1,nz+1},{0,0,4.5});
        connect("kmax",{1,1,nz+1,nx+1,ny+1,nz+1},{1,1,1,nx+1,ny+1,1},{0,0,-4.5});
        check(cg_close(file)); file=0;
        std::cout<<"hill cells="<<nx<<" x "<<ny<<" x "<<nz<<" = "<<static_cast<std::size_t>(nx)*ny*nz<<"; beta=2.2; x_cluster=0.45\n";
        return 0;
    } catch(const std::exception& e) { if(file) cg_close(file); std::cerr<<e.what()<<'\n'; return 1; }
}
