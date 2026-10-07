#include <cgnslib.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
void check(int x) { if(x) throw std::runtime_error(cg_get_error()); }
int main(int argc,char**argv) { try {
 if(argc!=3) return 2; int f,bases,zones,cd,pd,nb,nc; char name[33],donor[33]; cgsize_t size[9];
 check(cg_open(argv[1],CG_MODE_READ,&f));check(cg_nbases(f,&bases));check(cg_base_read(f,1,name,&cd,&pd));check(cg_nzones(f,1,&zones));check(cg_zone_read(f,1,1,name,size));
 if(bases!=1||zones!=1||cd!=3||pd!=3)throw std::runtime_error("dimensions");
 int nx=int(size[3]),ny=int(size[4]),nz=int(size[5]);const int ni=nx+1,nj=ny+1,nk=nz+1;
 std::cout<<std::setprecision(17)<<"CGNS read: zones="<<zones<<" nodes="<<ni<<"x"<<nj<<"x"<<nk<<" cells="<<nx<<"x"<<ny<<"x"<<nz<<" total="<<nx*ny*nz<<"\n";
 cgsize_t lo[3]={1,1,1},hi[3]={ni,nj,nk};std::vector<double>x(ni*nj*nk),y(x.size()),z(x.size());
 check(cg_coord_read(f,1,1,"CoordinateX",RealDouble,lo,hi,x.data()));check(cg_coord_read(f,1,1,"CoordinateY",RealDouble,lo,hi,y.data()));check(cg_coord_read(f,1,1,"CoordinateZ",RealDouble,lo,hi,z.data()));
 for(auto*v:{&x,&y,&z})std::cout<<"coordinate bounds="<<*std::min_element(v->begin(),v->end())<<":"<<*std::max_element(v->begin(),v->end())<<"\n";
 auto q=[&](int i,int j,int k){return(k*nj+j)*ni+i;};double minV=1e10,volume=0,mismatch=0;
 for(int k=0;k<nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){double dx=x[q(i+1,j,k)]-x[q(i,j,k)],dz=z[q(i,j,k+1)]-z[q(i,j,k)];double dy=.5*(y[q(i,j+1,k)]-y[q(i,j,k)]+y[q(i+1,j+1,k)]-y[q(i+1,j,k)]);double V=dx*dy*dz; if(!(V>0))throw std::runtime_error("volume");minV=std::min(minV,V);volume+=V;}
 for(int k=0;k<nk;++k)for(int j=0;j<nj;++j)mismatch=std::max(mismatch,std::abs(y[q(0,j,k)]-y[q(nx,j,k)]));
 for(int i=0;i<ni;++i)for(int j=0;j<nj;++j)mismatch=std::max(mismatch,std::abs(y[q(i,j,0)]-y[q(i,j,nz)]));
 std::cout<<"linear-cell min_volume="<<minV<<" total_volume="<<volume<<" periodic_coordinate_mismatch="<<mismatch<<"\n";
 check(cg_nbocos(f,1,1,&nb));check(cg_n1to1(f,1,1,&nc));if(nb!=2||nc!=4)throw std::runtime_error("boundary/connectivity counts");
 std::cout<<"physical_walls="<<nb<<" periodic_connections="<<nc<<"\n";
 for(int c=1;c<=nc;++c){cgsize_t r[6],d[6];int tr[3];float center[3],angle[3],shift[3];check(cg_1to1_read(f,1,1,c,name,donor,r,d,tr));check(cg_1to1_periodic_read(f,1,1,c,center,angle,shift));std::cout<<name<<" -> "<<donor<<" shift="<<shift[0]<<","<<shift[1]<<","<<shift[2]<<" transform="<<tr[0]<<","<<tr[1]<<","<<tr[2]<<"\n";}
 std::ofstream csv(argv[2]);csv<<std::setprecision(17)<<"i,j,x,y\n";for(int j=0;j<nj;++j)for(int i=0;i<ni;++i)csv<<i<<","<<j<<","<<x[q(i,j,0)]<<","<<y[q(i,j,0)]<<"\n";
 check(cg_close(f));return 0;
 }catch(const std::exception&e){std::cerr<<e.what();return 1;}}
