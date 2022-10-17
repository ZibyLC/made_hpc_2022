#include "matrix.h"

#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>

void Matrix::random_fill()
{

    for (int row = 0; row < rows(); ++row)
    {
        for (int col = 0; col < cols(); ++col)
        {
            std::random_device dev;
            std::mt19937 rng(dev());
            std::uniform_int_distribution<std::mt19937::result_type> dist(0,1000); // distribution in range [1, 6]

            mtrx[row][col] = dist(rng);
        }
    }
}

void Matrix::fill_row_with_same(int row, int k)
{

    for (int col = 0; col < mtrx.size(); ++col)
    {
        mtrx[row][col] = k;
    }
}

void Matrix::fill_col_with_same(int col, int k)
{

    for (int row = 0; row < mtrx.size(); ++row)
    {
        mtrx[row][col] = k;
    }
}

void Matrix::print() const
{
    std::cout << "Matrix of size (rows x cols) = " << rows() << "x" << cols() << ":" << std::endl;
    for (int row = 0; row < rows(); ++row)
    {
        for (int col = 0; col < cols(); ++col)
        {
            std::cout << std::setw(4) << mtrx[row][col] << ' ';
        }
        std::cout << std::endl;
    }
}

void print_column(const std::vector<int> &column)
{
    std::cout << "Column (vector) of size " << column.size() << ":" << std::endl;
    for (int i = 0; i < column.size(); ++i)
    {
        std::cout << std::setw(4) << column[i] << std::endl;
    }
    std::cout << std::endl;
}

unsigned long microsec()
{
    return std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

void test_square_matrix_x_square_matrix(int size)
{
    std::cout << "test_square_matrix_x_square_matrix() with size of:  " << size << std::endl;
    unsigned long start = microsec();

    Matrix first(size, size);

    first.random_fill();
    //first.print();

    Matrix second(size, size);

    second.random_fill();
    //second.print();

    Matrix result = first * second;

    unsigned long time_taken = microsec() - start;
    //result.print();
    std::cout << "Time taken: " << time_taken << " microseconds" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
}

void test_matrix_x_column()
{
    std::cout << "test_matrix_x_column()" << std::endl;
    unsigned long start = microsec();
    int rows = 3;
    int cols = 2;
    Matrix m(rows, cols);

    m.random_fill();
    m.print();

    std::vector<int> original_vector(cols, 1);
    print_column(original_vector);
    std::vector<int> result = m * original_vector;

    unsigned long time_taken = microsec() - start;
    print_column(result);
    std::cout << "Time taken: " << time_taken << " microseconds" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
}

void test_matrix_x_matrix_different_sizes()
{
    std::cout << "test_matrix_x_matrix_different_sizes()" << std::endl;
    unsigned long start = microsec();
    int rows = 3;
    int cols = 2;
    Matrix m(rows, cols);

    m.random_fill();
    m.print();

    std::vector<int> original_vector(cols, 1);
    print_column(original_vector);
    std::vector<int> result = m * original_vector;

    unsigned long time_taken = microsec() - start;
    print_column(result);
    std::cout << "Time taken: " << time_taken << " microseconds" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
}


int main()
{
    //test_matrix_x_column();
    //test_matrix_x_matrix_different_sizes();

    test_square_matrix_x_square_matrix(500);
    test_square_matrix_x_square_matrix(512);
    test_square_matrix_x_square_matrix(1000);
    test_square_matrix_x_square_matrix(1024);
    test_square_matrix_x_square_matrix(2000);
    test_square_matrix_x_square_matrix(2048);

    return 0;
}