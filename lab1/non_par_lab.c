#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

double f(double t, double x) {
    return x + t;
}

double phi(double x) {
    return cos(3.14159265 * x);
}

double psi(double t) {
    return exp(-t);
}

int idx(int k, int m, int M) {
    return k * M + m;
}

void fill_func(double(*f)(double, double), double* du_data, int M, int K, double tau, double h) {
for (int k = 0; k < K; k++) {
    for (int m = 0; m < M; m++) {
            int i = idx(k, m, M); 
            du_data[i] = f(k * tau, m * h);
        }
    }
}

void fill_init(double(*f)(double), double* arr, int N, double h) {
    for (int i = 0; i < N; ++i) {
        arr[i] = f(i * h);
    }
}

void file_print(const char* data_file_name, double* u, int M, int K) {
    FILE* data_file = fopen(data_file_name, "a");
    if (data_file) {
        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                int i = idx(k, m, M);
                fprintf(data_file, "%lf ", u[i]);
            }
            fprintf(data_file, "\n");
        }
        fclose(data_file);
    }
}

void krest(double* u, double* du_data, int M, int K, double tau, double h) {
    for (int m = 1; m < M; m++) {
        u[idx(1, m, M)] = du_data[idx(0, m, M)] * tau + (h - tau)/h * u[idx(0, m, M)] + tau / h * u[idx(0, m - 1, M)];
    }
    for (int k = 1; k < K; k++) {
        for (int m = 1; m < M; m++) {
            if (m < M - 1) {
                u[idx(k + 1, m, M)] = du_data[idx(k, m, M)] * 2 * tau + u[idx(k - 1, m, M)] + (u[idx(k, m - 1, M)] - u[idx(k, m + 1, M)]) * tau / h;
            }
            else {
                u[idx(k + 1, m, M)] = du_data[idx(k, m, M)] * tau + (h - tau) / h * u[idx(k, m, M)] + tau / h * u[idx(k, m - 1, M)];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    double T = 1;
    double X = 1;

    if (argc != 3) {
        fprintf(stderr, "[-] Usage %s K N\n", argv[0]);
        return -1;
    }
    int K = atoi(argv[1]);
    int M = atoi(argv[2]);
    
    double tau = T / K;
    double h = X / M;

    clock_t time_start = clock();

    double* du_data = (double*) calloc(K * M, sizeof(double));
    double* u     = (double*) calloc(K * M, sizeof(double));
    fill_func(f, du_data, M, K, tau, h);

    
    double* phi_data = (double*) calloc(K, sizeof(double));
    double* psi_data = (double*) calloc(M, sizeof(double));
    fill_init(psi, psi_data, K, tau); 
    fill_init(phi, phi_data, M, h); 
    
    for (int i = 0; i < M; ++i) {
        u[idx(0, i, M)] = phi_data[i];
    }
    
    for (int i = 0; i < K; ++i) {
        u[idx(i, 0, M)] = psi_data[i];
    }
    
    krest(u, du_data, M, K, tau, h);
    
    const char* data_file_name = "validate.txt";
    const char* file_name = "time.txt";
    file_print(data_file_name, u, M, K);

    clock_t time_end = clock();

    double total_time = (time_end - time_start) * 1000 / CLOCKS_PER_SEC;

    printf("Non parallel: %lf ms\n", total_time);
    FILE* data_file = fopen(file_name, "a");
    fprintf(data_file, "%lf ", total_time);
    fclose(data_file);

    free(psi_data);
    //free(phi_data);
    //free(du_data);
    //free(u);
    return 0;
}