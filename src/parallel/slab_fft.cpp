#include <wcns/parallel/slab_fft.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#if WCNS_HAS_FFTW
#include <fftw3.h>
#endif

namespace wcns {
namespace {
constexpr Real pi = 3.141592653589793238462643383279502884;
std::vector<int> offsets(const std::vector<int>& counts)
{
    std::vector<int> result(counts.size());
    long long total = 0;
    for (std::size_t r=0;r<counts.size();++r) {
        if (counts[r]<0 || total+counts[r]>std::numeric_limits<int>::max())
            throw MpiError("FFT MPI payload exceeds the int-count interface");
        result[r]=static_cast<int>(total); total+=counts[r];
    }
    return result;
}
}
std::vector<Real> exchange_peer_reals(const MpiRuntime& mpi, const std::vector<Real>& send,
                                     const std::vector<int>& counts,
                                     std::vector<int>* receive_counts)
{
    if (counts.size()!=static_cast<std::size_t>(mpi.size()))
        throw MpiError("FFT peer count size differs from communicator");
    auto sd=offsets(counts);
    if (static_cast<std::size_t>(sd.back()+counts.back())!=send.size())
        throw MpiError("FFT send payload size differs from peer counts");
    std::vector<int> rc=counts;
#if WCNS_HAS_MPI
    check_mpi(MPI_Alltoall(counts.data(),1,MPI_INT,rc.data(),1,MPI_INT,mpi.communicator()),
              "FFT Alltoall counts");
#endif
    auto rd=offsets(rc);
    std::vector<Real> receive(static_cast<std::size_t>(rd.back()+rc.back()));
#if WCNS_HAS_MPI
    Real dummy=0;
    check_mpi(MPI_Alltoallv(send.empty()?&dummy:send.data(),counts.data(),sd.data(),MPI_DOUBLE,
                            receive.empty()?&dummy:receive.data(),rc.data(),rd.data(),MPI_DOUBLE,
                            mpi.communicator()),"FFT Alltoallv payload");
#else
    receive=send;
#endif
    if(receive_counts) *receive_counts=std::move(rc);
    return receive;
}

struct SlabFft::Plans {
#if WCNS_HAS_FFTW
    fftw_complex* data=nullptr;
    fftw_plan forward=nullptr, inverse=nullptr;
    ~Plans() { if(forward) fftw_destroy_plan(forward); if(inverse) fftw_destroy_plan(inverse);
               if(data) fftw_free(data); }
#endif
};
SlabFft::SlabFft(const MpiRuntime& mpi,int n):mpi_(mpi),n_(n),plans_(std::make_unique<Plans>())
{
    // The radix-2 fallback and production FFTW path share the same supported contract.
    if(n<4 || n>512 || (n&(n-1))) throw std::invalid_argument("HIT FFT N must be a power of two in [4,512]");
#if WCNS_HAS_FFTW
    plans_->data=fftw_alloc_complex(n);
    if(!plans_->data) throw std::bad_alloc();
    plans_->forward=fftw_plan_dft_1d(n,plans_->data,plans_->data,FFTW_FORWARD,FFTW_ESTIMATE);
    plans_->inverse=fftw_plan_dft_1d(n,plans_->data,plans_->data,FFTW_BACKWARD,FFTW_ESTIMATE);
    if(!plans_->forward || !plans_->inverse) throw std::runtime_error("FFTW plan creation failed");
#endif
}
SlabFft::~SlabFft()=default;
const char* SlabFft::backend() const
{
#if WCNS_HAS_FFTW
    return "distributed-slab-MPI+FFTW3";
#else
    return "distributed-slab-MPI+radix2";
#endif
}
int SlabFft::owner(int x) const
{
    if(x<0 || x>=n_) throw std::out_of_range("FFT global x index");
    return static_cast<int>((static_cast<long long>(x+1)*mpi_.size()-1)/n_);
}
std::size_t SlabFft::index(int x,int y,int z) const
{ return (static_cast<std::size_t>(x)*n_+y)*n_+z; }
std::array<int,3> SlabFft::wave(std::size_t q) const
{
    std::array<int,3> k{{static_cast<int>(q/(n_*n_))+begin(),static_cast<int>(q/n_%n_),static_cast<int>(q%n_)}};
    for(auto& a:k) if(a>n_/2) a-=n_;
    return k;
}
void SlabFft::line(FourierField& a,bool inverse) const
{
#if WCNS_HAS_FFTW
    for(int i=0;i<n_;++i) { plans_->data[i][0]=a[i].real(); plans_->data[i][1]=a[i].imag(); }
    fftw_execute(inverse?plans_->inverse:plans_->forward);
    for(int i=0;i<n_;++i) a[i]={plans_->data[i][0],plans_->data[i][1]};
#else
    for(int i=1,j=0;i<n_;++i) {
        int bit=n_>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit;
        if(i<j) std::swap(a[i],a[j]);
    }
    for(int len=2;len<=n_;len*=2) {
        const FourierValue root=std::polar(Real{1},(inverse?2:-2)*pi/len);
        for(int i=0;i<n_;i+=len) {
            FourierValue w=1;
            for(int j=0;j<len/2;++j) { auto u=a[i+j],v=a[i+j+len/2]*w;
                a[i+j]=u+v; a[i+j+len/2]=u-v; w*=root; }
        }
    }
#endif
}
FourierField SlabFft::transpose(const FourierField& a) const
{
    std::vector<int> counts(mpi_.size());
    std::vector<Real> send; send.reserve(2*a.size());
    for(int r=0;r<mpi_.size();++r) {
        counts[r]=2*width()*width(r)*n_;
        for(int x=0;x<width();++x) for(int y=begin(r);y<begin(r+1);++y)
            for(int z=0;z<n_;++z) {const auto v=a[index(x,y,z)]; send.push_back(v.real());send.push_back(v.imag());}
    }
    auto receive=exchange_peer_reals(mpi_,send,counts);
    FourierField result(size()); std::size_t p=0;
    for(int r=0;r<mpi_.size();++r) for(int x=begin(r);x<begin(r+1);++x)
        for(int y=0;y<width();++y) for(int z=0;z<n_;++z) {
            result[index(y,x,z)]={receive[p],receive[p+1]};p+=2;
        }
    return result;
}
void SlabFft::transform(FourierField& a,bool inverse) const
{
    if(a.size()!=size()) throw std::invalid_argument("FFT local field size");
    FourierField work(n_);
    for(int x=0;x<width();++x) for(int y=0;y<n_;++y) {
        for(int z=0;z<n_;++z) work[z]=a[index(x,y,z)];
        line(work,inverse); for(int z=0;z<n_;++z) a[index(x,y,z)]=work[z];
    }
    for(int x=0;x<width();++x) for(int z=0;z<n_;++z) {
        for(int y=0;y<n_;++y) work[y]=a[index(x,y,z)];
        line(work,inverse); for(int y=0;y<n_;++y) a[index(x,y,z)]=work[y];
    }
    a=transpose(a);
    for(int y=0;y<width();++y) for(int z=0;z<n_;++z) {
        for(int x=0;x<n_;++x) work[x]=a[index(y,x,z)];
        line(work,inverse); for(int x=0;x<n_;++x) a[index(y,x,z)]=work[x];
    }
    a=transpose(a);
    if(!inverse) { const Real norm=Real{1}/(static_cast<Real>(n_)*n_*n_);for(auto& v:a) v*=norm; }
}
void SlabFft::forward(FourierField& a) const {transform(a,false);}
void SlabFft::inverse(FourierField& a) const {transform(a,true);}
} // namespace wcns
