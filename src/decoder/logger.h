/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#ifndef DECODER_LOGGER_h
#define DECODER_LOGGER_h

#include <array>
#include <iosfwd>
#include <sstream>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class LOGGER
{
public:
    constexpr static int MAX_VERBOSITY{3};

    using stream_type = std::stringstream;
    using stream_array = std::array<stream_type, MAX_VERBOSITY>;

    /*
     * User can increment/decrement this to adjust tabbing
     * of input to `stream()`.
     * */
    int tab_level{0};
private:
    /*
     * This logs debug info for all syndromes. Each array index
     * corresponds to a different verbosity level.
     * */
    stream_array info_strm_array_{};

    /*
     * This logs only for syndromes resulting in a logical
     * error.
     * */
    stream_type error_strm_;
public:
    LOGGER() =default;
    LOGGER(LOGGER&&) =default;

    stream_type& info(int verbosity);
    stream_type& error();

    void dump_info(std::ostream&, int verbosity);
    void dump_error(std::ostream&);

    void reset();
private:
    stream_type& prefix(stream_type&);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder

#endif // DECODER_LOGGER_h
