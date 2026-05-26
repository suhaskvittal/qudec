/*
 *  author: Suhas Vittal
 *  date:   25 May 2026
 * */

#include "globals.h"

#include <iostream>

std::ostream&
operator<<(std::ostream& ostrm, _die)
{
    ostrm << "\n";
    std::terminate();
    return ostrm;
}
