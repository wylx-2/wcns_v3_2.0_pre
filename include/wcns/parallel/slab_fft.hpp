#pragma once
#include <wcns/parallel/mpi_runtime.hpp>
#include <complex>
#include <array>
#include <memory>

namespace wcns {
using FourierValue = std::complex<Real>;
using FourierField = std::vector<FourierValue>;

// Equal-sided Cartesian box, arbitrary slab lengths (including empty slabs).
// Forward coefficients are normalized by N^3; inverse is unnormalized.
// No rank owns a replicated global field. The two slab transposes use Alltoallv.
class SlabFft {
public:
    SlabFft(const MpiRuntime&, int n);
    ~SlabFft();
    SlabFft(const SlabFft&) = delete;
    SlabFft& operator=(const SlabFft&) = delete;
    int n() const { return n_; }
    int begin(int rank) const { return n_ * rank / mpi_.size(); }
    int width(int rank) const { return begin(rank + 1) - begin(rank); }
    int begin() const { return begin(mpi_.rank()); }
    int width() const { return width(mpi_.rank()); }
    int owner(int x) const;
    std::size_t size() const { return static_cast<std::size_t>(width()) * n_ * n_; }
    std::size_t index(int local_x, int y, int z) const;
    std::array<int,3> wave(std::size_t index) const;
    void forward(FourierField&) const;
    void inverse(FourierField&) const;
    const char* backend() const;
private:
    struct Plans;
    void transform(FourierField&, bool inverse) const;
    FourierField transpose(const FourierField&) const;
    void line(FourierField&, bool inverse) const;
    const MpiRuntime& mpi_;
    int n_;
    std::unique_ptr<Plans> plans_;
};

// Peer payloads are concatenated by rank; counts are in Real elements.
std::vector<Real> exchange_peer_reals(const MpiRuntime&, const std::vector<Real>&,
                                     const std::vector<int>& send_counts,
                                     std::vector<int>* receive_counts = nullptr);
} // namespace wcns
