#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) 
{ 
    MPI_Comm server; 
    int buf = 150;
    char port_name[MPI_MAX_PORT_NAME]; 
 
    MPI_Init(&argc, &argv);

    MPI_Lookup_name("name", MPI_INFO_NULL, port_name);
 
    if (MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &server) != MPI_SUCCESS) {
        MPI_Finalize();
        return 1;
    }
    
    printf("Sending number: %d\n", buf);

    MPI_Send(&buf, 1, MPI_INT, 0, 1, server); 
    MPI_Comm_disconnect(&server); 
    MPI_Finalize(); 

    return 0; 
}