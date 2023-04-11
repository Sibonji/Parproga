#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{   
    MPI_Init(&argc, &argv);

    MPI_Comm intercomm;

    if (MPI_Comm_spawn("server", MPI_ARGV_NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE) != MPI_SUCCESS) {
        printf("MASTER: couldn't spawn server!\n");

        MPI_Finalize();
        return 0;
    }
    if (MPI_Comm_spawn("client", MPI_ARGV_NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE) != MPI_SUCCESS) {
        printf("MASTER: couldn't spawn client!\n");

        MPI_Finalize();
        return 0;
    }

    MPI_Finalize();
    return 0;
}