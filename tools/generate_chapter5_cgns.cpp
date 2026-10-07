#include <cgnslib.h>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <iomanip>
namespace {
void check(int s) {if(s!=CG_OK) throw std::runtime_error(cg_get_error());}
std::string name(int z) {return "Zone"+std::to_string(z+1);}
}
int main(int argc,char**argv) {
 int f=0;try {
    if(argc!=3) throw std::invalid_argument("usage: wcns_generate_chapter5_cgns plane.bin output.cgns");
    if(std::ifstream(argv[2]).good()) throw std::invalid_argument("refusing to overwrite mesh");
    std::ifstream in(argv[1],std::ios::binary);in.exceptions(std::ios::failbit|std::ios::badbit);
    char magic[4];in.read(magic,4);if(std::string(magic,4)!="XY25") throw std::runtime_error("bad plane magic");
    std::array<std::uint32_t,5> h{};in.read(reinterpret_cast<char*>(h.data()),20);
    const int kind=h[0],nb=h[1],nx=h[2],ny=h[3],nz=h[4];double span=0;in.read(reinterpret_cast<char*>(&span),8);
    if((kind!=0 && kind!=1)||nb<2||nb>128||nx<6||ny<6||nz<6||nx>2048||ny>1024||nz>512||!(span>0)) throw std::runtime_error("bad dimensions");
    const std::size_t plane=static_cast<std::size_t>(nx+1)*(ny+1),n=plane*(nz+1);
    if(n>5000000) throw std::runtime_error("zone too large for bounded generator");
    check(cg_set_file_type(CG_FILE_ADF));check(cg_open(argv[2],CG_MODE_WRITE,&f));int base=0;
    check(cg_base_write(f,"Base",3,3,&base));
    double minjac=1e300,maxskew=0;std::vector<std::vector<double>> left,right;
    for(int b=0;b<nb;++b) {
        std::vector<double> xp(plane),yp(plane),x(n),y(n),z(n);
        in.read(reinterpret_cast<char*>(xp.data()),plane*8);in.read(reinterpret_cast<char*>(yp.data()),plane*8);
        auto q=[&](int i,int j){return static_cast<std::size_t>(j)*(nx+1)+i;};
        std::vector<double> L,R;
        for(int j=0;j<=ny;++j) {L.push_back(xp[q(0,j)]);L.push_back(yp[q(0,j)]);R.push_back(xp[q(nx,j)]);R.push_back(yp[q(nx,j)]);}
        left.push_back(L);right.push_back(R);
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) {
            // Bilinear quadrilateral Jacobian at all four corners: stronger than area alone.
            for(int a=0;a<2;++a) for(int c=0;c<2;++c) {
                const double dx=xp[q(i+1,j+c)]-xp[q(i,j+c)],dy=yp[q(i+1,j+c)]-yp[q(i,j+c)];
                const double ex=xp[q(i+a,j+1)]-xp[q(i+a,j)],ey=yp[q(i+a,j+1)]-yp[q(i+a,j)];
                const double jac=(dx*ey-dy*ex)*span/nz;
                if(!(jac>0 && std::isfinite(jac))) throw std::runtime_error("nonpositive cell Jacobian in zone "+name(b)+" i="+std::to_string(i)+" j="+std::to_string(j));
                minjac=std::min(minjac,jac);maxskew=std::max(maxskew,std::abs((dx*ex+dy*ey)/(std::hypot(dx,dy)*std::hypot(ex,ey))));
            }
        }
        for(int k=0;k<=nz;++k) for(std::size_t p=0;p<plane;++p) {auto t=k*plane+p;x[t]=xp[p];y[t]=yp[p];z[t]=span*k/nz;}
        cgsize_t size[9]={nx+1,ny+1,nz+1,nx,ny,nz,0,0,0};int zone=0,coord=0;
        check(cg_zone_write(f,base,name(b).c_str(),size,Structured,&zone));
        check(cg_coord_write(f,base,zone,RealDouble,"CoordinateX",x.data(),&coord));
        check(cg_coord_write(f,base,zone,RealDouble,"CoordinateY",y.data(),&coord));
        check(cg_coord_write(f,base,zone,RealDouble,"CoordinateZ",z.data(),&coord));
        auto bc=[&](const char* n,BCType_t type,std::array<cgsize_t,6> r) {int id=0;check(cg_boco_write(f,base,zone,n,type,PointRange,2,r.data(),&id));};
        bc("wall",BCWallViscous,{1,1,1,nx+1,1,nz+1});bc("farfield",BCFarfield,{1,ny+1,1,nx+1,ny+1,nz+1});
        if(kind==1 && b==0) bc("inlet",BCInflow,{1,1,1,1,ny+1,nz+1});
        if(kind==1 && b==nb-1) bc("outlet",BCOutflow,{nx+1,1,1,nx+1,ny+1,nz+1});
    }
    double mismatch=0;
    for(int b=0;b<nb;++b) {
        const int zone=b+1;
        auto conn=[&](const char* n,int donor,std::array<cgsize_t,6> r,std::array<cgsize_t,6>d,double shift=0) {
            int tr[3]={1,2,3},id=0;check(cg_1to1_write(f,base,zone,n,name(donor).c_str(),r.data(),d.data(),tr,&id));
            if(shift!=0) {float zero[3]={0,0,0},trans[3]={0,0,static_cast<float>(shift)};check(cg_1to1_periodic_write(f,base,zone,id,zero,zero,trans));}
        };
        if(kind==0 || b>0) conn("imin",(b+nb-1)%nb,{1,1,1,1,ny+1,nz+1},{nx+1,1,1,nx+1,ny+1,nz+1});
        if(kind==0 || b<nb-1) {
            const int next=(b+1)%nb;for(std::size_t p=0;p<right[b].size();++p) mismatch=std::max(mismatch,std::abs(right[b][p]-left[next][p]));
            conn("imax",next,{nx+1,1,1,nx+1,ny+1,nz+1},{1,1,1,1,ny+1,nz+1});
        }
        conn("kmin",b,{1,1,1,nx+1,ny+1,1},{1,1,nz+1,nx+1,ny+1,nz+1},span);
        conn("kmax",b,{1,1,nz+1,nx+1,ny+1,nz+1},{1,1,1,nx+1,ny+1,1},-span);
    }
    if(mismatch>1e-12) throw std::runtime_error("interface mismatch");
    check(cg_close(f));f=0;
    std::cout<<std::setprecision(17)<<"zones="<<nb<<" cells_per_zone="<<nx<<"x"<<ny<<"x"<<nz<<" total_cells="<<static_cast<std::size_t>(nb)*nx*ny*nz
       <<" min_corner_jacobian="<<minjac<<" maximum_abs_cos_grid_angle="<<maxskew<<" interface_mismatch="<<mismatch<<"\n";
    return 0;
 }catch(const std::exception& e) {if(f) cg_close(f);std::cerr<<e.what()<<'\n';return 1;}
}
