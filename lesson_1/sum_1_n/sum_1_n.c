#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

int main(int argc, char* argv[])
{
	clock_t start, end;
	double cpu_time_used;
	
	start = clock();

	int pid, np,
		elements_per_process,
		n_elements_count;

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	int n = 1000000000;

	if (np == 1) {
		printf("Can't start programm with only one core!\n");
		
		MPI_Finalize();
		
		return 0;
	}

	// master process
	if (pid == 0) {
		int index, i;
		elements_per_process = n / (np - 1);
		int elem_to_count = 0;

		if (np > 1) {
			for (i = 0; i < np - 1; i++) {
				index = i * elements_per_process + 1;

				if (i == np - 2) {
					elem_to_count = n + 1 - index;
				}
				else
					elem_to_count = elements_per_process;

				MPI_Send(&elem_to_count,
						1, MPI_INT, i + 1, 0,
						MPI_COMM_WORLD);
				MPI_Send(&index,
						1,
						MPI_INT, i + 1, 0,
						MPI_COMM_WORLD);

				//printf("Elem to count: %d, start point: %d\n", elem_to_count, index);
			}
		}

		// collects partial sums from other processes
		float tmp;
		float sum = 0;
		for (i = 1; i < np; i++) {
			MPI_Recv(&tmp, 1, MPI_FLOAT,
					MPI_ANY_SOURCE, 0,
					MPI_COMM_WORLD,
					&status);
			int sender = status.MPI_SOURCE;

			sum += tmp;
		}

		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

		printf("Sum of array is : %f\n", sum);

		printf("Time taken: %f\n", cpu_time_used);
	}
	//slave process
	else {
		MPI_Recv(&n_elements_count,
				1, MPI_INT, 0, 0,
				MPI_COMM_WORLD,
				&status);

		int cur_num;
		MPI_Recv(&cur_num, 1,
				MPI_INT, 0, 0,
				MPI_COMM_WORLD,
				&status);

		//printf("Received elem to count: %d, start point: %d, pid: %d\n", n_elements_count, cur_num, pid);

		float partial_sum = 0;
		for (int i = cur_num; i < (n_elements_count + cur_num); i++)
			partial_sum += (float)1 / i;

		//printf("Partial sum is: %f, pid: %d\n", partial_sum, pid);

		MPI_Send(&partial_sum, 1, MPI_FLOAT,
				0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;
}