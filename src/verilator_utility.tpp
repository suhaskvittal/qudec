/*
 *  author: Suhas Vittal
 *  date:   19 June 2026
 * */

#include <algorithm>
#include <cassert>
#include <numeric>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

constexpr size_t
verilator_width(size_t N, size_t W)
{
    constexpr size_t vl_wide_width{32};
    const size_t total_width = N*W;
    return (total_width+vl_wide_width-1) / vl_wide_width;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <size_t N, size_t W, class IterType, class Callback> VWideDataType<N,W>
verilator_build_packed_array(IterType begin, IterType end, const Callback& f_cb)
{
    using WordType = uint64_t;
    constexpr size_t word_bit_width = 8*sizeof(WordType);
    constexpr size_t K_WORD = (W+word_bit_width-1)/word_bit_width;
    constexpr size_t K_B32 = (sizeof(WordType)/sizeof(uint32_t))*K_WORD;

    VWideDataType<N,W> out{};

    size_t word_idx{0}, 
           bit_idx{0};
    for (auto it = begin; it != end; it++)
    {
        auto repr = f_cb(*it);
        // convert to 32b representation
        auto r32b = std::bit_cast<std::array<uint32_t, K_B32>>(repr);
        // update `out`
        size_t bits_applied_in_current_word{0},
               bits_applied{0};
        while (bits_applied < W)
        {
            size_t bits_to_apply = std::min(W - bits_applied, size_t{32} - bits_applied_in_current_word);
            bits_to_apply = std::min(bits_to_apply, size_t{32} - bit_idx);

            size_t i = bits_applied / 32;
            uint32_t mask = (bits_to_apply == 32) 
                                    ? std::numeric_limits<uint32_t>::max()
                                    : (1<<bits_to_apply)-1;

            uint32_t& x = r32b[i];
            out[word_idx] |= (x & mask) << bit_idx;
            x >>= bits_to_apply;
            bit_idx += bits_to_apply;
            if (bit_idx == 32)
            {
                word_idx++;
                bit_idx = 0;
            }
            bits_applied += bits_to_apply;
            bits_applied_in_current_word = (bits_applied_in_current_word + bits_to_apply) % 32;
        }
    }
    return out;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
