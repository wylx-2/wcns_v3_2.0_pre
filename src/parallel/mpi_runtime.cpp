#include <wcns/parallel/mpi_runtime.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <utility>

namespace wcns {

#if WCNS_HAS_MPI

void check_mpi(int status, const char* operation)
{
    if (status == MPI_SUCCESS) {
        return;
    }
    std::array<char, MPI_MAX_ERROR_STRING> message {};
    int length = 0;
    MPI_Error_string(status, message.data(), &length);
    throw MpiError(std::string(operation) + ": "
                   + std::string(message.data(), static_cast<std::size_t>(length)));
}

#endif

MpiRuntime::MpiRuntime(int& argc, char**& argv)
{
#if WCNS_HAS_MPI
    int finalized = 0;
    check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
    if (finalized != 0) {
        throw MpiError("MPI has already been finalized");
    }

    int initialized = 0;
    check_mpi(MPI_Initialized(&initialized), "MPI_Initialized");
    if (initialized == 0) {
        check_mpi(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level_),
                  "MPI_Init_thread");
        owns_mpi_ = true;
    } else {
        check_mpi(MPI_Query_thread(&thread_level_), "MPI_Query_thread");
    }
    if (thread_level_ < MPI_THREAD_FUNNELED) {
        if (owns_mpi_) {
            MPI_Finalize();
            owns_mpi_ = false;
        }
        throw MpiError("MPI does not provide MPI_THREAD_FUNNELED");
    }

    int raw_rank = 0;
    check_mpi(MPI_Comm_rank(MPI_COMM_WORLD, &raw_rank), "MPI_Comm_rank");
    check_mpi(MPI_Comm_size(MPI_COMM_WORLD, &size_), "MPI_Comm_size");
    check_mpi(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN),
              "MPI_Comm_set_errhandler");
    rank_ = static_cast<RankId>(raw_rank);
    enabled_ = true;
#else
    static_cast<void>(argc);
    static_cast<void>(argv);
#endif
}

MpiRuntime::~MpiRuntime()
{
#if WCNS_HAS_MPI
    if (!owns_mpi_) {
        return;
    }
    int finalized = 0;
    if (MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0) {
        MPI_Finalize();
    }
#endif
}

void MpiRuntime::barrier() const
{
#if WCNS_HAS_MPI
    check_mpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier");
#endif
}

Real MpiRuntime::sum(Real local_value) const
{
#if WCNS_HAS_MPI
    Real global_value = 0.0;
    check_mpi(MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD),
              "MPI_Allreduce sum");
    return global_value;
#else
    return local_value;
#endif
}

Real MpiRuntime::min(Real local_value) const
{
#if WCNS_HAS_MPI
    Real global_value = 0.0;
    check_mpi(MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD),
              "MPI_Allreduce min");
    return global_value;
#else
    return local_value;
#endif
}

Real MpiRuntime::max(Real local_value) const
{
#if WCNS_HAS_MPI
    Real global_value = 0.0;
    check_mpi(MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD),
              "MPI_Allreduce max");
    return global_value;
#else
    return local_value;
#endif
}

bool MpiRuntime::all_true(bool local_value) const
{
#if WCNS_HAS_MPI
    const int local = local_value ? 1 : 0;
    int global = 0;
    check_mpi(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD),
              "MPI_Allreduce logical and");
    return global != 0;
#else
    return local_value;
#endif
}

bool MpiRuntime::all_equal(std::uint64_t local_value) const
{
#if WCNS_HAS_MPI
    unsigned long long local = static_cast<unsigned long long>(local_value);
    unsigned long long minimum = 0;
    unsigned long long maximum = 0;
    check_mpi(MPI_Allreduce(&local, &minimum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD),
              "MPI_Allreduce uint64 minimum");
    check_mpi(MPI_Allreduce(&local, &maximum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD),
              "MPI_Allreduce uint64 maximum");
    return minimum == maximum;
#else
    static_cast<void>(local_value);
    return true;
#endif
}

std::string MpiRuntime::broadcast_string(std::string value, RankId root) const
{
    if (root < 0 || root >= size_) {
        throw MpiError("broadcast root is outside MPI rank range");
    }
#if WCNS_HAS_MPI
    unsigned long long size = rank_ == root ? static_cast<unsigned long long>(value.size()) : 0;
    check_mpi(MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, static_cast<int>(root), MPI_COMM_WORLD),
              "MPI_Bcast string size");
    if (size > static_cast<unsigned long long>(std::numeric_limits<int>::max())) {
        throw MpiError("broadcast string exceeds MPI int count range");
    }
    if (rank_ != root) value.resize(static_cast<std::size_t>(size));
    if (size > 0) {
        check_mpi(MPI_Bcast(value.data(),
                            static_cast<int>(size),
                            MPI_CHAR,
                            static_cast<int>(root),
                            MPI_COMM_WORLD),
                  "MPI_Bcast string data");
    }
#else
    static_cast<void>(root);
#endif
    return value;
}

std::vector<Real> MpiRuntime::gather_reals(std::vector<Real> local_values, RankId root) const
{
    if (root < 0 || root >= size_) {
        throw MpiError("gather root is outside MPI rank range");
    }
    if (local_values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw MpiError("local real gather count exceeds MPI int range");
    }
#if WCNS_HAS_MPI
    const int local_count = static_cast<int>(local_values.size());
    std::vector<int> counts;
    if (rank_ == root) counts.resize(static_cast<std::size_t>(size_));
    check_mpi(MPI_Gather(&local_count,
                         1,
                         MPI_INT,
                         rank_ == root ? counts.data() : nullptr,
                         1,
                         MPI_INT,
                         static_cast<int>(root),
                         MPI_COMM_WORLD),
              "MPI_Gather real counts");
    std::vector<int> displacements;
    std::vector<Real> result;
    int valid_total = 1;
    if (rank_ == root) {
        displacements.resize(static_cast<std::size_t>(size_));
        std::size_t total = 0;
        for (int rank = 0; rank < size_; ++rank) {
            displacements[static_cast<std::size_t>(rank)] = static_cast<int>(total);
            total += static_cast<std::size_t>(counts[static_cast<std::size_t>(rank)]);
            if (total > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                valid_total = 0;
                break;
            }
        }
        if (valid_total != 0) {
            result = std::move(local_values);
            result.resize(total);
            const auto root_offset
                = static_cast<std::size_t>(displacements[static_cast<std::size_t>(root)]);
            if (root_offset != 0 && local_count != 0) {
                std::move_backward(result.begin(),
                                   result.begin() + local_count,
                                   result.begin()
                                       + static_cast<std::ptrdiff_t>(
                                           root_offset + static_cast<std::size_t>(local_count)));
            }
        }
    }
    check_mpi(MPI_Bcast(&valid_total, 1, MPI_INT, static_cast<int>(root), MPI_COMM_WORLD),
              "MPI_Bcast real gather validity");
    if (valid_total == 0) {
        throw MpiError("global real gather count exceeds MPI int range");
    }
    check_mpi(MPI_Gatherv(rank_ == root ? MPI_IN_PLACE
                                        : (local_values.empty() ? nullptr : local_values.data()),
                          rank_ == root ? 0 : local_count,
                          MPI_DOUBLE,
                          rank_ == root && !result.empty() ? result.data() : nullptr,
                          rank_ == root ? counts.data() : nullptr,
                          rank_ == root ? displacements.data() : nullptr,
                          MPI_DOUBLE,
                          static_cast<int>(root),
                          MPI_COMM_WORLD),
              "MPI_Gatherv real values");
    return result;
#else
    static_cast<void>(root);
    return local_values;
#endif
}

std::vector<Real> MpiRuntime::all_gather_reals(std::vector<Real> local_values) const
{
    if (local_values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw MpiError("local real all-gather count exceeds MPI int range");
    }
#if WCNS_HAS_MPI
    const int local_count = static_cast<int>(local_values.size());
    std::vector<int> counts(static_cast<std::size_t>(size_));
    check_mpi(MPI_Allgather(&local_count,
                            1,
                            MPI_INT,
                            counts.data(),
                            1,
                            MPI_INT,
                            MPI_COMM_WORLD),
              "MPI_Allgather real counts");
    std::vector<int> displacements(static_cast<std::size_t>(size_));
    std::size_t total = 0;
    for (int rank = 0; rank < size_; ++rank) {
        displacements[static_cast<std::size_t>(rank)] = static_cast<int>(total);
        const auto count = static_cast<std::size_t>(counts[static_cast<std::size_t>(rank)]);
        if (total > static_cast<std::size_t>(std::numeric_limits<int>::max()) - count) {
            throw MpiError("global real all-gather count exceeds MPI int range");
        }
        total += count;
    }
    std::vector<Real> result(total);
    check_mpi(MPI_Allgatherv(local_values.empty() ? nullptr : local_values.data(),
                             local_count,
                             MPI_DOUBLE,
                             result.empty() ? nullptr : result.data(),
                             counts.data(),
                             displacements.data(),
                             MPI_DOUBLE,
                             MPI_COMM_WORLD),
              "MPI_Allgatherv real values");
    return result;
#else
    return local_values;
#endif
}

std::vector<Real> MpiRuntime::scatter_reals(const std::vector<Real>& root_values,
                                            const std::vector<std::size_t>& counts,
                                            RankId root) const
{
    if (root < 0 || root >= size_) {
        throw MpiError("scatter root is outside MPI rank range");
    }
    std::vector<int> integer_counts;
    std::vector<int> displacements;
    int valid = 1;
    if (rank_ == root) {
        if (counts.size() != static_cast<std::size_t>(size_)) valid = 0;
        integer_counts.resize(static_cast<std::size_t>(size_));
        displacements.resize(static_cast<std::size_t>(size_));
        std::size_t total = 0;
        if (valid != 0) {
            for (int rank = 0; rank < size_; ++rank) {
                const auto count = counts[static_cast<std::size_t>(rank)];
                if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())
                    || total > static_cast<std::size_t>(std::numeric_limits<int>::max()) - count) {
                    valid = 0;
                    break;
                }
                integer_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
                displacements[static_cast<std::size_t>(rank)] = static_cast<int>(total);
                total += count;
            }
        }
        if (valid != 0 && total != root_values.size()) valid = 0;
    }
#if WCNS_HAS_MPI
    check_mpi(MPI_Bcast(&valid, 1, MPI_INT, static_cast<int>(root), MPI_COMM_WORLD),
              "MPI_Bcast real scatter validity");
    if (valid == 0) throw MpiError("real scatter counts are invalid");
    int local_count = 0;
    check_mpi(MPI_Scatter(rank_ == root ? integer_counts.data() : nullptr,
                          1,
                          MPI_INT,
                          &local_count,
                          1,
                          MPI_INT,
                          static_cast<int>(root),
                          MPI_COMM_WORLD),
              "MPI_Scatter real counts");
    std::vector<Real> result(static_cast<std::size_t>(local_count));
    check_mpi(MPI_Scatterv(rank_ == root && !root_values.empty() ? root_values.data() : nullptr,
                           rank_ == root ? integer_counts.data() : nullptr,
                           rank_ == root ? displacements.data() : nullptr,
                           MPI_DOUBLE,
                           result.empty() ? nullptr : result.data(),
                           local_count,
                           MPI_DOUBLE,
                           static_cast<int>(root),
                           MPI_COMM_WORLD),
              "MPI_Scatterv real values");
    return result;
#else
    if (valid == 0) throw MpiError("real scatter counts are invalid");
    return root_values;
#endif
}

} // namespace wcns
