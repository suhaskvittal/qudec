#ifndef DECODER_TYPES_h
#define DECODER_TYPES_h

#include <stim.h>

namespace decoder
{

using SyndromeType = stim::simd_bits<64>;
using SyndromeRef = stim::simd_bits_range_ref<64>;
using ObsType = SyndromeType;
using ObsRef = SyndromeRef;

} // namespace decoder

#endif
