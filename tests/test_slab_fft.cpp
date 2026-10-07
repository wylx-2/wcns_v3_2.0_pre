#include <wcns/parallel/slab_fft.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv){
 wcns::MpiRuntime mpi(argc,argv);
 try {
  for(int n:{8,16}) {
   wcns::SlabFft fft(mpi,n);wcns::FourierField a(fft.size());double physical=0;
   const double pi=std::acos(-1.);
   for(int x=0;x<fft.width();++x)for(int y=0;y<n;++y)for(int z=0;z<n;++z){
    double v=1.7+std::sin(2*pi*(fft.begin()+x+2*y+3*z)/n)+.25*std::cos(4*pi*z/n);
    a[fft.index(x,y,z)]=v;physical+=v*v;}
   auto original=a;fft.forward(a);double spectral=0,error=0;
   for(std::size_t q=0;q<a.size();++q){spectral+=std::norm(a[q]);auto k=fft.wave(q);
    if(k[0]==0 && k[1]==0 && k[2]==0)error=std::abs(a[q]-wcns::FourierValue(1.7,0));}
   spectral=mpi.sum(spectral);physical=mpi.sum(physical)/(double(n)*n*n);
   if(std::abs(spectral-physical)>1e-12 || mpi.max(error)>1e-12)throw std::runtime_error("FFT Parseval or zero mode failed");
   fft.inverse(a);error=0;for(std::size_t q=0;q<a.size();++q)error=std::max(error,std::abs(a[q]-original[q]));
   if(mpi.max(error)>1e-12)throw std::runtime_error("FFT inverse/transpose failed");
   if(mpi.rank()==0)std::cout<<fft.backend()<<" N="<<n<<" ranks="<<mpi.size()<<" Parseval="<<spectral-physical<<" inverse="<<error<<'\n';
  }
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
