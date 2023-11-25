#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    int N = 0;
    if (argc == 2) {
        N = atoi(argv[1]);
    }
    else {
        printf("Usage:\n./sum_1_N [N]\nWhere N >= 1!\n");
        return 1;
    }

    if (N <= 0) {
        printf("N need to be more than zero!\n");
        return 1;
    }

    float sum = 0;
    float partial_sum = 0;

    #pragma omp parallel private(partial_sum) shared(sum)
    {
        #pragma omp for
            for(int i = 1; i <= N; i++){
                partial_sum += 1/(float)i;
            }
        
        #pragma omp critical
        {
            sum += partial_sum;
        }
    }

    printf("Sum is: %f\n", sum);

    sum = 0;
    for (float i = 1; i <= N; i++) {
        sum += 1 / i;
    }

    printf("Expected is: %f\n", sum);

    return 0;
}