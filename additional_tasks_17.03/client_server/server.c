#include "mpi.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) 
{ 
    MPI_Comm client; 
    MPI_Status status; 
    char port_name[MPI_MAX_PORT_NAME]; 
    int buf;
    int size;
 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    if (size != 1) {
        printf("Server is too big!\n");
        MPI_Finalize();
        return 0;
    }
    if (MPI_Open_port(MPI_INFO_NULL, port_name) != MPI_SUCCESS){
        MPI_Finalize();
        return 0;
    }
    MPI_Publish_name("name", MPI_INFO_NULL, port_name);

    printf("server available at %s\n", port_name); 
    while (1) { 
        MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &client); 
        MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status); 
            
        printf("Server received number: %d!\n", buf);

        MPI_Comm_free(&client); 
            
        MPI_Close_port(port_name); 

        MPI_Unpublish_name("name", MPI_INFO_NULL, port_name);
        MPI_Finalize(); 
        return 0;
    } 
}