/* Author: OpenAI Codex. Shared entry point for qudec and qudec_chains. */
#include "qudec_cli.h"
#include <iostream>
#include <exception>
#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

int main(int argc, char* argv[])
{
    int world_rank = 0;
#if defined(ENABLE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
    int status = 0;
    try {
        status = run_cli(argc, argv, world_rank,
#if defined(QUDEC_CHAINS)
            true
#else
            false
#endif
        );
    } catch (const std::exception& error) {
        if (world_rank == 0) std::cerr << argv[0] << ": " << error.what() << '\n';
        status = 1;
#if defined(ENABLE_MPI)
        MPI_Abort(MPI_COMM_WORLD, status);
#endif
    }
#if defined(ENABLE_MPI)
    MPI_Finalize();
#endif
    return status;
}

