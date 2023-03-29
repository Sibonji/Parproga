#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    int pid, np;

    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int msg = 0;


    MPI_Send(&msg, 1, MPI_INT, (pid + 1) % np, 0, MPI_COMM_WORLD);
    printf("Sent message to %d, message is: %d, my id is: %d\n", (pid + 1) % np, msg, pid);

    int res;

    if (pid != 0) {
        MPI_Recv(&res, 1, MPI_INT, pid - 1, 0, MPI_COMM_WORLD, &status);
        printf("Received message from %d, message is: %d, my id is: %d\n", pid - 1, res, pid);
    }
    else {
        MPI_Recv(&res, 1, MPI_INT, np - 1, 0, MPI_COMM_WORLD, &status);
        printf("Received message from %d, message is: %d, my id is: %d\n", np - 1, res, pid);
    }

    MPI_Finalize();

    return 0;
}