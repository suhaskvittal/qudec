/*
 *  author: Suhas Vittal
 *  date:   19 June 2026
 * */

#ifndef VERILATOR_UTILITY_h
#define VERILATOR_UTILITY_h

#include <verilated.h>

#include <cstdint>
#include <cstddef>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
 * Returns the width of `VLWide` required for `N` values with width `W`
 * */
constexpr size_t verilator_width(size_t N, size_t W);

template <size_t N, size_t W> 
using VWideDataType = VlWide<verilator_width(N,W)>;

/*
 * Builds a packed array from an iterator range and
 * a given callback. The callback is called on every
 * element in the iterator range and should return a
 * `std::array<uint64_t, K>` that corresponds to
 * that entry's binary representation. `K` is equal
 * to the number of 64-bit words required to represent
 * the entry. For example, if W <= 64, then K = 1.
 * If 64 < W <= 128, then K = 2. Et cetera.
 * */
template <size_t N, size_t W, class IterType, class Callback>
VWideDataType<N,W> verilator_build_packed_array(IterType begin, IterType end, const Callback&);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#include "verilator_utility.tpp"

#endif
