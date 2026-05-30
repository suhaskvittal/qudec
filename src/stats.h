/*
 *  author: Suhas Vittal
 *  date:   25 May 2026
 * */

#ifndef STATS_h
#define STATS_h

#include <limits>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class T>
class STATS_HISTOGRAM
{
public:
    const T range_min;
    const T range_max;
    const T bucket_width;
    const size_t num_buckets;
private:
    T min_{ std::numeric_limits<T>::max() };
    T max_{};
    T sum_{};
    T sum_of_sq_{};
    
    std::vector<size_t> buckets_;
    size_t count_{0};
public:
    STATS_HISTOGRAM(T range_min, T range_max, size_t num_buckets);

    template <class U> void add(U);

    /*
     * Getting stats:
     * */
    double mean() const;
    double std() const;
    T min() const { return min_; }
    T max() const { return max_; }

    size_t operator[](size_t idx) const { return buckets_[idx]; };
    size_t underflow_count() const { return buckets_[underflow_idx()]; }
    size_t overflow_count() const { return buckets_[overflow_idx()]; }

    std::string to_string_some() const;
    std::string to_string_full() const;
private:
    size_t underflow_idx() const { return num_buckets+1; }
    size_t overflow_idx() const { return num_buckets+2; }
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <class T>
void print_stat(std::ostream&, std::string_view name, T);

template <class T>
std::ostream& operator<<(std::ostream&, const STATS_HISTOGRAM<T>&);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#include "stats.tpp"

#endif // STATS_h
