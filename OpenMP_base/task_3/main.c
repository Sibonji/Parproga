#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv){
    int i;

    #pragma omp parallel
    {
        for( int i = 0; i < omp_get_max_threads(); i++){
            if(i == omp_get_thread_num()){
                printf("Hello from process: %d\n", omp_get_thread_num());
            }
            #pragma omp barrier
        }
    }

    return 0;
}
