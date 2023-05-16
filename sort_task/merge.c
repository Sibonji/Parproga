#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#define OUTPUT_NUM 10
#define MASTER 0

void merge(int *, int *, int, int, int);
void mergeSort(int *, int *, int, int);

int main(int argc, char** argv) {
	int process_rank;
	int num_processes;
    double timer_start;
    double timer_end;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	int array_size = atoi(argv[1]) * pow(2, 20) / sizeof(int);
	int *array = calloc(array_size, sizeof(int));
	
	int c;
	srand(time(NULL));
	for(c = 0; c < array_size; c++) {	
		array[c] = rand() % array_size;
    }

    if (process_rank == MASTER) {
        printf("\n%d numbers in start array:\n",OUTPUT_NUM);
        for (int j = 0; j < array_size; j++) {
            if ((j % (array_size / OUTPUT_NUM)) == 0) {
                printf("%d ",array[j]);
            }
        } 
        printf("\n\n");

        timer_start = MPI_Wtime();
    }

	int size = array_size / num_processes;

	int *sub_array = calloc(size, sizeof(int));
    MPI_Scatter(array, size, MPI_INT, sub_array, size, MPI_INT, 0, MPI_COMM_WORLD);
	
	int *tmp_array = calloc(size, sizeof(int));
	mergeSort(sub_array, tmp_array, 0, (size - 1));
	
	int *sorted = NULL;
	if(process_rank == MASTER) {
		sorted = calloc(array_size, sizeof(int));	
	}
	
	MPI_Gather(sub_array, size, MPI_INT, sorted, size, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(process_rank == MASTER) {
		int *other_array = calloc(array_size, sizeof(int));
		mergeSort(sorted, other_array, 0, (array_size - 1));
		
        printf("\n%d numbers in sorted array:\n",OUTPUT_NUM);
        for (int j = 0; j < array_size; j++) {
            if ((j % (array_size / OUTPUT_NUM)) == 0) {
                printf("%d ",sorted[j]);
            }
        } 
        printf("\n\n");

        timer_end = MPI_Wtime();

        printf("Time Elapsed (Sec): %f\n", timer_end - timer_start);

		free(sorted);
		free(other_array);		
	}

	free(array);
	free(sub_array);
	free(tmp_array);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

    return 0;
}

void merge(int *a, int *b, int l, int m, int r) {
	int h, i, j, k;
	h = l;
	i = l;
	j = m + 1;
	
	while((h <= m) && (j <= r)) {
		if(a[h] <= a[j]) {
			b[i] = a[h];
			h++;	
	    }	
		else {
			b[i] = a[j];
			j++;	
		}

		i++;	
	}
		
	if(m < h) {
		for(k = j; k <= r; k++) {
			b[i] = a[k];
			i++;	
		}		
	}	
	else {
		for(k = h; k <= m; k++) {
			b[i] = a[k];
			i++;
		}	
	}
		
	for(k = l; k <= r; k++) {
		a[k] = b[k];
	}

    return;	
}

void mergeSort(int *a, int *b, int l, int r) {
	int m;
	
	if(l < r) {
		m = (l + r)/2;
		
		mergeSort(a, b, l, m);
		mergeSort(a, b, (m + 1), r);
		merge(a, b, l, m, r);		
	}

    return;
}