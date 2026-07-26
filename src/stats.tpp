/*
 *  author: Suhas Vittal
 *  date: 25 May 2026
 * */

#define TEMPL_PARAM template <class T>
#define TEMPL_CLASS StatsHistogram<T>

#include <iomanip>
#include <iostream>
#include <sstream>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM
TEMPL_CLASS::StatsHistogram(T _range_min,
                             T _range_max,
                             size_t _num_buckets)
    :range_min(_range_min),
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

TEMPL_PARAM double
TEMPL_CLASS::mean() const
{
    return static_cast<double>(sum_) / static_cast<double>(count_);
}

TEMPL_PARAM double
TEMPL_CLASS::std() const
{
    double mean_of_sq = static_cast<double>(sum_of_sq_) / static_cast<double>(count_);
    return std::sqrt(mean_of_sq - mean());
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

TEMPL_PARAM std::string
TEMPL_CLASS::to_string_some() const
{
    std::stringstream strm;
    strm << "mean=" << mean() << ", std=" << std() << ", min=" << min() << ", max=" << max();
    return strm.str();
}

TEMPL_PARAM std::string
TEMPL_CLASS::to_string_full() const
{
    return "";
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
