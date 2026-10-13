/*
 *  author: Suhas Vittal
 *  date:   13 October 2026
 * */

#ifndef DECODER_SC_PDZ_h
#define DECODER_SC_PDZ_h

/*
 * `PDZ` = Predecoderization. This is essentially the concept of
 * taking the outputs of an initial decode and using them to seed
 * decisions of another decoder.
 *
 * For surface codes, this can be done by taking the matchings produced
 * by a decoder and eliminating the "obvious" matchings to reduce the 
 * amount of work done by another decoder.
 * */

namespace dec
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

result_type pdz_apply_to_syndrome(SyndromeRef, const result_type& first_dec_result);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace dec

#endif
