#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include "mpi.h"

int main(int argc, char* argv[]) {
    int N = atoi(argv[1]);

    int size = 0;
    int rank = 0;

    double value = 3.141592;

    int res = MPI_Init(&argc, &argv);
    if (res != MPI_SUCCESS) {
        perror("[-] Error starting MPI program. Programm was terminated.\n");
        MPI_Abort(MPI_COMM_WORLD, res);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 1) {
        double start = MPI_Wtime();
        double recv = 0;
        for (int i = 0; i < N; ++i) {
            MPI_Send(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&recv, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        double end = MPI_Wtime();

        double total_time = (end - start) * 1000; // ms
        printf("Single send time: %f ms\n", total_time / 2 / N);
        printf("Total time: %f ms\n", total_time);
    }

    if (rank == 0) {
        double recv = 0;
        for (int i = 0; i < N; ++i) {
            MPI_Recv(&recv, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&recv, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            assert(recv == value);
        }

    }

    MPI_Finalize();

    return 0;
}