#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{   
    MPI_Init(&argc, &argv);

    char port_name[MPI_MAX_PORT_NAME];
    MPI_Lookup_name("name", MPI_INFO_NULL, port_name);

    MPI_Comm server;
    if (MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &server) != MPI_SUCCESS) {
        printf("CLIENT: couldn't connect to the port!\n");
        
        MPI_Finalize();
        return 0;
    }

    int num = 3;
    printf("CLIENT: sending number %d\n", num);
    MPI_Send(&num, 1, MPI_INT, 0, 0, server);

    MPI_Finalize();
    return 0;
}