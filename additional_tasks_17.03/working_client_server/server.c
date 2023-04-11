#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{  
    MPI_Init(&argc, &argv);

    // Open port.
    char port_name[MPI_MAX_PORT_NAME];
    if (MPI_Open_port(MPI_INFO_NULL, port_name) != MPI_SUCCESS) {
        printf("SERVER: couldn't open port!\n");
        MPI_Finalize();
        return 0;
    }

    printf("SERVER: opened port with name: %s\n", port_name);

    MPI_Comm client;
    MPI_Publish_name("name", MPI_INFO_NULL, port_name);
    MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &client);

    int num;
    MPI_Recv(&num, 1, MPI_INT, 0, 0, client, MPI_STATUS_IGNORE);
    printf("SERVER: received number = %d\n", num);

    MPI_Unpublish_name("name", MPI_INFO_NULL, port_name);
    MPI_Finalize();

    return 0;
}