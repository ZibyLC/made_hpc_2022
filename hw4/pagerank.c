#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <omp.h>
#include <math.h>

const double EPS = 0.0001;

void norm_cols_matrix(const size_t N, double *matrix)
{
    double sum[N];
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        sum[i] = 0;

#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for (size_t j = 0; j < N; ++j)
                sum[i] += matrix[j * N + i];
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for (size_t j = 0; j < N; ++j)
            if (abs(sum[j]) > EPS)
                matrix[i * N + j] /= sum[j];
    return;
}


double norm2(const size_t N, const double *vector)
{
    double norm = 0;
#pragma omp parallel for reduction(+: norm)
    for (size_t i = 0; i < N; ++i)   
        norm += vector[i] * vector[i];
    return sqrt(norm);
}


double norm2_vector(const size_t N, double *vector)
{
    const double norm = norm2(N, &vector[0]);

    if (abs(norm) > EPS)
    {
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        vector[i] /= norm;
    }
    return norm;
}


void mult_matrix_on_vector(const size_t N, const double *MATRIX, const double *VECTOR, double *res)
{
    double temp[N];
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        temp[i] = 0;
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        for(size_t k = 0; k < N; ++k)
            temp[i] += MATRIX[i * N + k] * VECTOR[k];
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        res[i] = temp[i];
    return;
}


void pagerank(const size_t N, double *matrix, double *x)
{
    // initial approximation
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)   
        x[i] = 1;

    double norm_cur, norm_last;
    norm_last = norm2_vector(N, &x[0]);

    printf("pagerank iters \n norm = %f \n", norm_last);

    // norm graph columns
    norm_cols_matrix(N, matrix);
    
    // 1st iteration
    mult_matrix_on_vector(N, &matrix[0], &x[0], &x[0]);
    norm_cur = norm2_vector(N, &x[0]);

    // iteration process
    for (; abs(norm_cur - norm_last) >= EPS;)  
    {
        printf("norm = %f\n", norm_cur);
        norm_last = norm_cur;
        mult_matrix_on_vector(N, &matrix[0], &x[0], &x[0]);
        norm_cur = norm2_vector(N, &x[0]);
    }   
    printf("norm = %f\n", norm_cur);

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


void print_vector(const size_t N, const double *vector)
{
    for (size_t i = 0; i < N; ++i)
        printf("%f ", vector[i]); 
    printf("\n"); 
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
    const size_t N = (argc > 1) ? atoi(argv[1]) : 5;
	double matrix[N * N];
	double ranks[N];
    double naiveranks[N];
    const size_t MAX_LINK_COUNT = (argc > 2) ? atoi(argv[2]) : 1;

#pragma omp parallel private(tid, seed)
    {
        tid = omp_get_thread_num();
        seed = seeds[tid];
        srand(seed);
// generate rnd graph
#pragma omp for 
        for (size_t i = 0; i < N; ++i)
        {    
            for (size_t j = 0; j < N; ++j)
                if (i != j)
                    matrix[i * N + j] = rand() % (MAX_LINK_COUNT + 1);
                else
                    matrix[i * N + j] = 0;
        }
    }

// print generated graph
    printf("\n graph \n");
    print_matrix(N, matrix);
    printf("\n");

//  calculate naive ranks
#pragma omp parallel for 
    for (size_t i = 0; i < N; ++i)
    {    
        naiveranks[i] = 0;
        for (size_t j = 0; j < N; ++j)
            naiveranks[i] += matrix[i * N + j];      
    }

// pagerank algo
    pagerank(N, &matrix[0], &ranks[0]);
    printf("\n pagerank & ranks \n");
    print_vector(N, &ranks[0]);

// naive ranking
    printf("\n naive ranking & inner links cnt \n");
    print_vector(N, &naiveranks[0]);
    printf("\n");

    return 0;
}
