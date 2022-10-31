#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <omp.h>


void mult_matrix(const size_t N, const double *MATRIX_1, const double *MATRIX_2, double *res)
{
    double temp[N * N];
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for (size_t j = 0; j < N; ++j)
            temp[i * N + j] = 0;
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for (size_t j = 0; j < N; ++j)
            for(size_t k = 0; k < N; ++k)
				temp[i * N + j] += MATRIX_1[i * N + k] * MATRIX_2[k * N + j];
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for (size_t j = 0; j < N; ++j)
            res[i * N + j] = temp[i * N + j];
    return;
}


void mpow(const size_t k, size_t *p2, size_t *k2) 
{
    p2[0] = 0;
    k2[0] = 1;

    if (k < 2)
        return;

    while (k2[0] <= k) 
    {
        k2[0] <<= 1;
        ++p2[0];
    }

    if (k2[0] > k)
    {
        k2[0] >>= 1;
        --p2[0];
    }
    return;
}


void pow_matrix(const size_t k, const size_t N, double *matrix)
{
    if (k == 0)
    {
#pragma omp parallel for 
        for (size_t i = 0; i < N; ++i)   
            for (size_t j = 0; j < N; ++j)
                if (i != j)
                    matrix[i * N + j] = 0;
                else
                    matrix[i * N + j] = 1;
        return;
    }
    size_t k2 = 0, p2 = 0;
    double temp[N * N];
    mpow(k, &p2, &k2);
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for (size_t j = 0; j < N; ++j)
            temp[i * N + j] = matrix[i * N + j];

    for (size_t i = 0; i < p2; ++i)
        mult_matrix(N, &matrix[0], &matrix[0], &matrix[0]);

    for (size_t i = 0; i < (k - k2); ++i)
        mult_matrix(N, &temp[0], &matrix[0], &matrix[0]);
    return;
}


void print_matrix(const size_t N, const double *matrix)
{
    for (size_t i = 0; i < N; ++i)
    {    
        for (size_t j = 0; j < N; ++j)
            printf("%f ", matrix[i * N + j]); 
        printf("\n"); 
    }
    return;
}


void seedThreads(const size_t nThreads, unsigned int* seeds) {
    unsigned int seed;
	int thread_id;
    #pragma omp parallel private (seed, thread_id)
    {
		unsigned int seed = (unsigned) time(NULL);
        thread_id = omp_get_thread_num();
        seeds[thread_id] = (seed & 0xFFFFFFF0) | (thread_id + 1);
    }
    return;
}


int main(int argc, char* argv[])
{
    const size_t nThreads = (argc > 3) ? atoi(argv[3]) : 4;

    srand(time(NULL));
    unsigned int seeds[nThreads];
    omp_set_num_threads(nThreads);
    seedThreads(nThreads, seeds);

    size_t tid, seed;
    const size_t N = (argc > 1) ? atoi(argv[1]) : 3;
    const size_t POW = (argc > 2) ? atoi(argv[2]) : 3;
	double matrix[N * N];
    const int MAX_VALUE = 10;

#pragma omp parallel private(tid, seed)
    {
        tid = omp_get_thread_num();
        seed = seeds[tid];
        srand(seed);
// generate rnd matrix
#pragma omp for 
        for (size_t i = 0; i < N; ++i)
        {    
            for (size_t j = 0; j < N; ++j)
            { 
                matrix[i * N + j] = rand() % MAX_VALUE;
            }
        }
    }

// print generated matrix
    printf("\n matrix \n");
    print_matrix(N, matrix);

// pow matrix
    pow_matrix(POW, N, &matrix[0]);
    
// print matrix
    printf("\n matrix ^ %d \n", POW);
    print_matrix(N, matrix);
    return 0;
}