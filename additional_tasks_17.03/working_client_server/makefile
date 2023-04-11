all:
	mpicc master.c -o master_program
	mpicc server.c -o server
	mpicc client.c -o client
	mpirun -n 1 master_program