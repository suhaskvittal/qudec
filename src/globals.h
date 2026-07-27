/*
 *  author: Suhas Vittal
 *  date:   25 May 2026
 * */

#ifndef GLOBALS_h
#define GLOBALS_h

#include <iosfwd>

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

struct _die{};
std::ostream& operator<<(std::ostream&, _die);

/*
 * `mpi_datatype<T>()` maps a C++ arithmetic type to its `MPI_Datatype`. Selecting on
 * size/signedness (rather than exact type identity) avoids the `size_t` vs `uint64_t`
 * mismatch that bites on platforms where they are distinct types (e.g. macOS:
 * `unsigned long` vs `unsigned long long`).
 *
 * In a non-MPI build there is no MPI type system, so calling this is a programming error:
 * the function reports to `std::cerr` and terminates.
 * */
#if defined(ENABLE_MPI)
template <class T> MPI_Datatype mpi_datatype();
#else
template <class T> void mpi_datatype();
#endif

#include "globals.tpp"

#endif // GLOBALS_h
