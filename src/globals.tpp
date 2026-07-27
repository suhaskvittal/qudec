/*
 *  author: Suhas Vittal
 *  date:   25 May 2026
 * */

#include <iostream>
#include <type_traits>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#if defined(ENABLE_MPI)

template <class T> MPI_Datatype
mpi_datatype()
{
    if constexpr (std::is_same_v<T, float>)
        return MPI_FLOAT;
    else if constexpr (std::is_same_v<T, double>)
        return MPI_DOUBLE;
    else if constexpr (std::is_integral_v<T> && std::is_signed_v<T>)
    {
        if constexpr (sizeof(T) == 8) return MPI_INT64_T;
        else if constexpr (sizeof(T) == 4) return MPI_INT32_T;
        else if constexpr (sizeof(T) == 2) return MPI_INT16_T;
        else return MPI_INT8_T;
    }
    else
    {
        static_assert(std::is_integral_v<T>, "unsupported type for mpi_datatype()");
        if constexpr (sizeof(T) == 8) return MPI_UINT64_T;
        else if constexpr (sizeof(T) == 4) return MPI_UINT32_T;
        else if constexpr (sizeof(T) == 2) return MPI_UINT16_T;
        else return MPI_UINT8_T;
    }
}

#else

template <class T> void
mpi_datatype()
{
    std::cerr << "mpi_datatype() called in a non-MPI build" << _die{};
}

#endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
