#include <thread>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <mpi.h>

double a = 0;
double b = 1.99;

double f(double x){
	return cos(1 / (2 - x));
}

double integrate (double leftBorder, double rightBorder, double precision) {
    double dx = sqrt(precision);
    double sum = 0;
    double prevsum = sum;
    double i;
    
    do {
        prevsum = sum;
        sum = 0;
        i = leftBorder;
    
        while (i < rightBorder)
        {
            sum = sum + dx * f(i + dx);
            i += dx;
        }
        
        dx /= 2;
        //printf("%.10lf\n", fabs(prevsum - sum));
    } while (fabs(prevsum - sum) > precision);

    return sum;
}

int main(int argc, char* argv[]) 
{
    double epsilon = 0;
    if(argc == 2){
        epsilon = atof(argv[1]);
    }
    else {
        epsilon = 0.00001;
	}

    int my_rank;
	int p;

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

    double startWTime;
	if(my_rank == 0)
	{
		startWTime = MPI_Wtime();
	}

	double local_a = a + my_rank * (b-a) / p;
	double local_b = a + (my_rank + 1) * (b-a) / p;

    double sum = 0;
    sum = integrate(local_a, local_b, epsilon);

    double total;
    MPI_Reduce(&sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0 , MPI_COMM_WORLD);

    if(my_rank ==0){
		double endWTime = MPI_Wtime();
		std::cout<<"Integral " << total << std::endl;
		std::cout<<"Number of processors used = "<<p<<std::endl;
		std::cout<<"Time elapsed: "<<(endWTime-startWTime)*1000<<"ms"<<std::endl;
	}

	MPI_Finalize();
    
    return 0;
}