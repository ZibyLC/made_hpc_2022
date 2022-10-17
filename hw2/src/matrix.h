#include <vector>


class Matrix
{
public:
    Matrix(int rows, int cols);

    void random_fill(); // matrix will be filled with random values for testing
    void fill_row_with_same(int row, int k); // fill all elements of a specified row with the same value
    void fill_col_with_same(int col, int k); // fill all elements of a specified column with the same value

    std::vector<int> operator *(const std::vector<int> &column) const;
    Matrix operator *(const Matrix &rhs) const;

    void print() const;

    int rows() const;
    int cols() const;

private:
    std::vector<std::vector<int>> mtrx;

};