#ifndef DECODER_TYPES_h
#define DECODER_TYPES_h

#include <stim.h>

namespace decoder
{

using syndrome_type = stim::simd_bits<64>;
using syndrome_ref = stim::simd_bits_range_ref<64>;
using obs_type = syndrome_type;
using obs_ref = syndrome_ref;

} // namespace decoder

#endif
