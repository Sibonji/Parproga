#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

//Функция ДУ
double du(double t, double x) {
    return x + t;
}

// Функция, задающая начальное условие на границе t = 0
double phi(double x) {
    return cos(3.1415926535 * x);
}

// Функция, задающая начальное условие на границе x = 0
double psi(double t) {
    return exp(-t);
}

// Функция возвращает
// k -- первый индекс двумерного массива (номер строки)
// m -- второй индекс двумерного массива (номер столбца)
// M -- количество столбцов в двумерном массиве
// a[k][m]
int GetIdx(int k, int m, int M) {
    return k * M + m;
}

// Функция заполняет значения функции в узлах сетки
// f -- указатель на функцию ДУ
// du_data -- указатель на массив значений функии в узлах сетки
// M -- количество точек на пространственном диапазоне
// K -- количество точек на временном диапазоне
// tau -- шаг сетки временного диапазона
// h -- шаг сетки пространственного диапазона
void FillFunctionValues(double(*func)(double, double), double* du_data, int M, int K, double tau, double h) {
for (int k = 0; k < K; ++k) {
    for (int m = 0; m < M; ++m) {
            int idx = GetIdx(k, m, M); // индекс одномерного массива
            du_data[idx] = func(k * tau, m * h);
        }
    }
}

// Функция заполняет узлы сетки на границе в соответствии с начальными условиями
// f -- указатель на функцию, задающую начальное значение на границе
// arr -- указатель на массив значений на границе
// N -- количество точек соответствующего диапазона (временного или пространственного)
// h -- шаг сетки 
void FillInitialValues(double(*func)(double), double* arr, int N, double h) {
    for (int i = 0; i < N; ++i) {
        arr[i] = func(i*h);
    }
}

void PutData2File(const char* data_file_name, double* u, int M, int K) {
    FILE* fd = fopen(data_file_name, "a");
    if (fd) {
        for (int m = 0; m < M; ++m) {
            for (int k = 0; k < K; ++k) {
                int idx = GetIdx(k, m, M);
                fprintf(fd, "%lf ", u[idx]);
            }
            fprintf(fd, "\n");
        }
        fclose(fd);
    }
}

int GetRankStart(int rank, int size, int N) {
    int n_range = N / size;
    int mod_size = N % size;

    int rank_start = 0;
    if (rank < mod_size) {
        rank_start = rank * (n_range + 1);
    } 
    else {
        rank_start = rank*n_range + mod_size;
    }

    return rank_start;
}

int GetRankEnd(int rank, int size, int N) {
    int n_range = N / size;
    int mod_size = N % size;

    int rank_end = 0;
    int rank_start = GetRankStart(rank, size, N);
    if (rank < mod_size) {
        rank_end = rank_start + n_range;
    } 
    else {
        rank_end = rank_start + n_range - 1;
    }

    return rank_end;
}

void PutData2FileParallel(int rank, int size, const char* data_file_name, double* u, int M, int K) {
    if (rank != 0) {
        int special_signal = 0;
        MPI_Recv(&special_signal, 1, MPI_INT, rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    PutData2File(data_file_name, u, M, K);
    if (rank != size-1) {
        int special_signal = 69;
        MPI_Send(&special_signal, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }
}

void CrossScheme(int rank, int size, double* u, double* f_arr, int rank_M, int K, int M, double tau, double h) {
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
    // обсчитываем первый слой схемой уголка (k = 0)
    //------------------------------------------------------------------------------------------
    for (int m = 1; m < rank_M; ++m)
        u[GetIdx(1, m, rank_M)] = 
                                    f_arr[GetIdx(0, m + rank*rank_M, M)] * tau + 
                                    (h - tau) / h * u[GetIdx(0, m, rank_M)] + 
                                    tau / h * u[GetIdx(0, m - 1, rank_M)];  

    // каждый процесс, кроме последнего отправляет своему соседу крайнюю
    // правую точку своего промежутка
    if (rank != size - 1)
        MPI_Isend(&u[GetIdx(0, rank_M - 1, rank_M)], 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &request);

    // Каждый процесс, кроме нулевого должен получить от предыдущего процесса 
    // крайнюю левую точку для старта расчета своего промежутка
    if (rank != 0) {
        double u_prev = 0.0;
        MPI_Irecv(&u_prev, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &request); //destination process receives
        MPI_Wait(&request, &status);
        u[GetIdx(1, 0, rank_M)] = 
                                    f_arr[GetIdx(0, 0 + rank * rank_M, M)] * tau + 
                                    (h - tau) / h * u[GetIdx(0, 0, rank_M)] + 
                                    tau / h * u_prev;
    }
    //------------------------------------------------------------------------------------------

    for (int k = 1; k < K-1; ++k) {
        for (int m = 1; m < rank_M; ++m) {
            if (m < rank_M-1) 
                u[GetIdx(k + 1, m, rank_M)] = 
                                        f_arr[GetIdx(k, m + rank*rank_M, M)] * 2 * tau + 
                                        u[GetIdx(k-1, m, rank_M)] + 
                                        (
                                            u[GetIdx(k, m-1, rank_M)] - u[GetIdx(k, m+1, rank_M)]
                                        ) * tau / h;
            else // точки на левой границе считаем схемой левого уголка
                u[GetIdx(k + 1, m, rank_M)] = 
                                        f_arr[GetIdx(k, m + rank * rank_M, M)] * tau + 
                                        (h - tau) / h * u[GetIdx(k, m, rank_M)] + 
                                        tau / h * u[GetIdx(k, m-1, rank_M)];
        }

        if (rank != size - 1) {
            MPI_Isend(&u[GetIdx(k, rank_M - 1, rank_M)], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request);
        }

        if (rank != 0) {
            double u_prev = 0.0;
            MPI_Irecv(&u_prev, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &request); //destination process receives
            MPI_Wait(&request, &status);
            u[GetIdx(k + 1, 0, rank_M)] = 
                                        f_arr[GetIdx(k, 0 + rank * rank_M, M)] * 2 * tau + 
                                        u[GetIdx(k - 1, 0, rank_M)] + 
                                        (
                                            u_prev - u[GetIdx(k, 1, rank_M)]
                                        ) * tau / h;
        }
    }
}

int main(int argc, char* argv[]) {
    /*
    Задаем T, X, K, M -- 
    границы пространственного и временного диапазона 
    и количество точек на расчетной сетке.
    */
    double T = 1;
    double X = 1;

    /*Ввод количества точек на расчетной сетке через параметры командной строки*/
    if (argc != 3) {
        fprintf(stderr, "[-] Usage mpiexec -n 4 %s K N\n", argv[0]);
        return -1;
    }
    //int K = atoi(argv[1]);
    //int M = atoi(argv[2]);
    int K = 1000;
    int M = 1000;

    // шаг сетки по времени и по пространству
    double tau = T / K;
    double h = X / M;

    /*
    Двумерные массивы arr[K][M] будем хранить в виде одномерных массивов
    arr[K*M]. В этом случае будет удобнее передавать их в функции по указателям и,
    в случае надобности, выделять динамически.
    */
    double* du_data = (double*) calloc(K*M, sizeof(double));
    //Заполнить массив значений функции в узлах сетки
    FillFunctionValues(du, du_data, M, K, tau, h);

    // начальные условия задачи
    double* phi_data = (double*) calloc(K, sizeof(double));
    double* psi_data = (double*) calloc(M, sizeof(double));
    7
    //Задаем начальные условия
    FillInitialValues(psi, psi_data, K, tau); // значения функции psi(t) в узлах сетки
    FillInitialValues(phi, phi_data, M, h); // значения функции phi(t) в узлах сетки

    int size = 0;
    int rank = 0;
    // double max_runk_time = 0;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double time_start = MPI_Wtime();

    // устанавливаем начало и конец диапазона, который обрабатывает текущий процесс
    int rank_M_start = GetRankStart(rank, size, M);
    int rank_M_end = GetRankEnd(rank, size, M);
    int rank_M = abs(rank_M_end - rank_M_start + 1);

    double* u = (double*) calloc(K * rank_M, sizeof(double));
    for (int i = rank_M_start; i <= rank_M_end; ++i) {
        u[GetIdx(0, i-rank_M*rank, rank_M)] = phi_data[i];
    }

    if (rank == 0) {
        for (int i = 0; i < K; ++i)
            u[GetIdx(i, 0, rank_M)] = psi_data[i];
    }

    CrossScheme(rank, size, u, du_data, rank_M, K, M, tau, h); // Крест

    const char* data_file_name = "data.txt";
    PutData2FileParallel(rank, size, data_file_name, u, rank_M, K);
    
    double time_end = MPI_Wtime();
    double rank_time = (time_end - time_start) * 1000;
    
    if (rank == size - 1) {
        const char* time_file_name = "time.txt";
        FILE* time_out_file = fopen(time_file_name, "a");
        if (time_out_file) 
            fprintf(time_out_file, "%lf ", rank_time);
        fclose(time_out_file);
    }

    MPI_Finalize();

    free(u);
    free(psi_data);
    free(phi_data);
    free(du_data);
    return 0;
}