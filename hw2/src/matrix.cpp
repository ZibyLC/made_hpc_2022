#include "matrix.h"

#include <stdexcept>

Matrix::Matrix(int rows, int cols)
{
    mtrx.resize(rows);

    for (int row = 0; row < rows; ++row)
    {
        mtrx[row].resize(cols);
    }
}

int Matrix::rows() const
{
    return mtrx.size();
}

int Matrix::cols() const
{
    return mtrx[0].size();
}

std::vector<int> Matrix::operator *(const std::vector<int> &column) const
{
    if (column.size() != cols())
        throw std::invalid_argument("amount of columns in the matrix should be equal to amount of elements in the column!");

    std::vector<int> result(rows());

    // elements of each line of the matrix should be multiplied by each element of the column
    for (int row = 0; row < rows(); ++row)
    {
        for (int col = 0; col < cols(); ++col)
        {
            result[row] += mtrx[row][col] * column[col];
        }
    }
    return result;
}

Matrix Matrix::operator *(const Matrix &other) const
{
    if (other.rows() != cols())
        throw std::invalid_argument("amount of columns in first matrix should be equal to amount of rows in second matrix!");

    Matrix result(rows(), other.cols());

	for(int i = 0; i < result.rows(); ++i)
    {
		for(int j = 0; j < result.cols(); ++j)
        {
			for(int k = 0; k < cols(); ++k)
            {
				result.mtrx[i][j] += mtrx[i][k] * other.mtrx[k][j];
            }
        }
    }

    return result;
}