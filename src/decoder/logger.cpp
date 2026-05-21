/*
 *  author: Suhas Vittal
 *  date:   16 May 2026
 * */

#include "decoder/logger.h"

#include <cassert>
#include <iostream>

namespace decoder
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace 
{

using stream_type = LOGGER::stream_type;

} // anon

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stream_type&
LOGGER::info(int v)
{
    assert(v < MAX_VERBOSITY);
    auto& strm = info_strm_array_[v];
    return prefix(strm);
}

stream_type&
LOGGER::error()
{
    return prefix(error_strm_);
}

void
LOGGER::dump_info(std::ostream& os, int v)
{
    os << info_strm_array_[v].str();
}

void
LOGGER::dump_error(std::ostream& os)
{
    os << error_strm_.str();
}

void
LOGGER::reset()
{
    for (auto& strm : info_strm_array_)
        strm.str("");
    error_strm_.str("");
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stream_type&
LOGGER::prefix(stream_type& strm)
{
    for (int i = 0; i < tab_level; i++)
        strm << "\t";
    return strm;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

} // namespace decoder
