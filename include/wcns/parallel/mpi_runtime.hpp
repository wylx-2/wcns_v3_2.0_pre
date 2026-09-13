#pragma once

#include <wcns/core/types.hpp>

#include <stdexcept>
#include <cstdint>
#include <string>
#include <vector>

#if WCNS_HAS_MPI
#include <mpi.h>
#endif

namespace wcns {

class MpiError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MpiRuntime {
public:
    MpiRuntime(int& argc, char**& argv);
    ~MpiRuntime();

    MpiRuntime(const MpiRuntime&) = delete;
    MpiRuntime& operator=(const MpiRuntime&) = delete;
    MpiRuntime(MpiRuntime&&) = delete;
    MpiRuntime& operator=(MpiRuntime&&) = delete;

    [[nodiscard]] bool enabled() const noexcept { return enabled_; }
    [[nodiscard]] RankId rank() const noexcept { return rank_; }
    [[nodiscard]] int size() const noexcept { return size_; }
    [[nodiscard]] int thread_level() const noexcept { return thread_level_; }

    void barrier() const;
    [[nodiscard]] Real sum(Real local_value) const;
    [[nodiscard]] Real min(Real local_value) const;
    [[nodiscard]] Real max(Real local_value) const;
    [[nodiscard]] bool all_true(bool local_value) const;
    [[nodiscard]] bool all_equal(std::uint64_t local_value) const;
    [[nodiscard]] std::string broadcast_string(std::string value, RankId root = 0) const;
    [[nodiscard]] std::vector<Real> gather_reals(std::vector<Real> local_values,
                                                 RankId root = 0) const;
    [[nodiscard]] std::vector<Real> all_gather_reals(std::vector<Real> local_values) const;
    [[nodiscard]] std::vector<Real> scatter_reals(const std::vector<Real>& root_values,
                                                  const std::vector<std::size_t>& counts,
                                                  RankId root = 0) const;

#if WCNS_HAS_MPI
    [[nodiscard]] MPI_Comm communicator() const noexcept { return MPI_COMM_WORLD; }
#endif

private:
    bool enabled_ = false;
    bool owns_mpi_ = false;
    RankId rank_ = 0;
    int size_ = 1;
    int thread_level_ = 0;
};

#if WCNS_HAS_MPI
void check_mpi(int status, const char* operation);
#endif

} // namespace wcns
