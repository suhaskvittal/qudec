/*
 *  author: Suhas Vittal
 *  date: 25 May 2026
 * */

#define TEMPL_PARAM template <class T>
#define TEMPL_CLASS StatsHistogram<T>

#include "globals.h"

#include <iomanip>
#include <iostream>
#include <sstream>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM
TEMPL_CLASS::StatsHistogram(std::string_view _name, 
                             T _range_min,
                             T _range_max,
                             size_t _num_buckets)
    :name(_name),
    range_min(_range_min),
    range_max(_range_max),
    bucket_width((range_max-range_min) / _num_buckets),
    num_buckets(_num_buckets),
    buckets_(num_buckets+2, 0)
{
    // verify that `bucket_width*num_buckets = (range_max-range_min)`
    if (bucket_width * num_buckets != (range_max-range_min))
        std::cerr << "STATS_HISTOGRAM: bucket_width*num_buckets != range_max-range_min" << _die{};
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM template <class U> void
TEMPL_CLASS::add(U _x)
{
    T x = static_cast<T>(_x);

    size_t idx;
    if (x < range_min)
        idx = underflow_idx();
    else if (x >= range_max)
        idx = overflow_idx();
    else
        idx = (x - range_min) / bucket_width;
    buckets_[idx]++;
    min_ = std::min(x, min_);
    max_ = std::max(x, max_);
    sum_ += x;
    sum_of_sq_ += x*x;
    count_++;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM void
TEMPL_CLASS::mpi_accumulate()
{
#if defined(ENABLE_MPI)
    const MPI_Datatype t_dtype = mpi_datatype<T>();
    const MPI_Datatype count_dtype = mpi_datatype<size_t>();

    // Sum per-bucket counts (including the underflow/overflow slots) across all ranks.
    MPI_Allreduce(MPI_IN_PLACE, buckets_.data(), static_cast<int>(buckets_.size()),
                    count_dtype, MPI_SUM, MPI_COMM_WORLD);

    // Reduce the running accumulators: min/max are global extrema, the rest are sums.
    MPI_Allreduce(MPI_IN_PLACE, &min_,       1, t_dtype,     MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &max_,       1, t_dtype,     MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_,       1, t_dtype,     MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_of_sq_, 1, t_dtype,     MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &count_,     1, count_dtype, MPI_SUM, MPI_COMM_WORLD);
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM double
TEMPL_CLASS::mean() const
{
    return static_cast<double>(sum_) / static_cast<double>(count_);
}

TEMPL_PARAM double
TEMPL_CLASS::std() const
{
    double mean_of_sq = static_cast<double>(sum_of_sq_) / static_cast<double>(count_);
    return std::sqrt(mean_of_sq - mean()*mean());
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM std::string
TEMPL_CLASS::to_string_some() const
{
    std::stringstream strm;
    strm << name << " : mean=" << mean() << ", std=" << std() << ", min=" << min() << ", max=" << max();
    return strm.str();
}

TEMPL_PARAM std::string
TEMPL_CLASS::to_string_full() const
{
    std::stringstream strm;
    strm << name << " =================================\n"
            << "mean = " << mean() 
            << ", std = " << std() 
            << ", min = " << min() 
            << ", max = " << max() 
            << ", count = " << count_
            << "\n";

    strm << "<" << range_min << ":\t" << underflow_count() << "\n";
    for (size_t i = 0; i < num_buckets; i++)
    {
        T from = range_min + bucket_width*i,
          to = range_min + bucket_width*(i+1);
        strm << from << " <= X < " << to << ":\t" << buckets_[i] << "\n";
    }
    strm << ">=" << range_max << ":\t" << overflow_count() << "\n";
    return strm.str();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class T> void
print_stat(std::ostream& ostrm, std::string_view name, T val)
{
    ostrm << std::setw(48) << std::left << name;
    if constexpr (std::is_floating_point<T>::value)
        ostrm << std::setprecision(5);
    ostrm << std::setw(16) << std::right << val << "\n";
}

template <class T> std::ostream&
operator<<(std::ostream& ostrm, const StatsHistogram<T>& hist)
{
    return (ostrm << hist.to_string_some());
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#undef TEMPL_PARAM
#undef TEMPL_CLASS
