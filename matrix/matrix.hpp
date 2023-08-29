#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

template <typename T>
class Matrix
{
    int rows{}, cols{};
    T **data = nullptr;

public:
    Matrix(const Matrix &rhs);
    Matrix(Matrix &&rhs) noexcept;
    Matrix(){};

    Matrix &operator=(const Matrix &rhs);
    Matrix &operator=(Matrix &&rhs) noexcept;

    ~Matrix();

public:
    Matrix(int cols, int rows, T val = {});
    Matrix(std::initializer_list<std::initializer_list<T>> list);

    static Matrix eye(int n, int m);

public: // operations
    Matrix &negate() &;
    Matrix &transpose() &;

public:
    int ncols() const;
    int nrows() const;

    T trace() const;

    // The det method using the sum of the products
    // of the elements of any one row or column and their cofactors.
    T det() const;

private:
    T det(Matrix<T> matrix) const;

public:
    void swap(Matrix &rhs) noexcept;

    bool operator==(const Matrix<T> &rhs) const;

    // Compares matrix1 and matrix2 ((rows1 == rows2) && (cols1 == cols2)),
    // element by element,
    // returning rows x cols matrix containing 1 where predicate is true, and 0 otherwise.
    Matrix equal(const Matrix &other) const;

    // Compares matrix1 and matrix2 ((rows1 == rows2) && (cols1 == cols2)),
    // element by element,
    // returning rows x cols matrix containing 1 where predicate is true, and 0 otherwise.
    Matrix less(const Matrix &other) const;

    void dump(std::ostream &os) const;

private:
    struct ProxyRow
    {
        T *row;
        const T &operator[](int n) const { return row[n]; }
        T &operator[](int n) { return row[n]; }

    public:
        ProxyRow(T *h) : row(h) {}
    };

public:
    ProxyRow operator[](int m) { return ProxyRow(data[m]); }
    const ProxyRow operator[](int m) const { return ProxyRow(data[m]); }
};

template <typename T>
int Matrix<T>::ncols() const { return cols; }

template <typename T>
int Matrix<T>::nrows() const { return rows; }

template <typename T>
Matrix<T> Matrix<T>::eye(int n, int m)
{
    assert(n == m);

    Matrix<T> result(n, m, 0);
 
    for (int i = 0; i < m; ++i)
        result[i][i] = 1;

    return result;
}

template <typename T>
T Matrix<T>::det() const { assert(cols == rows); return det(*this); }

template <typename T>
T Matrix<T>::det(Matrix<T> matrix) const
{
    int sz_ = matrix.ncols();

    if (sz_ == 1)
        return matrix[0][0];

    else if (sz_ == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    else
    {
        T resultD = 0;
        for (int k = 0; k < sz_; k++)
        {
            Matrix<T> subMatrix(sz_ - 1, sz_ - 1);

            for (int i = 1; i < sz_; i++)
            {
                for (int j = 0; j < sz_; j++)
                {
                    if (j == k)
                        continue;
                    else if (j < k)
                        subMatrix[i - 1][j] = matrix[i][j];
                    else
                        subMatrix[i - 1][j - 1] = matrix[i][j];
                }
            }
            resultD += std::pow(-1, k + 2) * matrix[0][k] * det(subMatrix);
        }
        return resultD;
    }
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T> &rhs) const
{
    assert((cols == rhs.ncols()) && (rows == rhs.nrows()));
    
    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            if (data[i][k] != rhs[i][k])
                return false;

    return true;
}

template <typename T>
Matrix<T> Matrix<T>::less(const Matrix &other) const
{
    assert((cols == other.ncols()) && (rows == other.nrows()));

    Matrix<T> result(cols, rows, 0);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            if (data[i][k] < other[i][k])
                result[i][k] = 1;

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::equal(const Matrix &other) const
{
    assert((cols == other.ncols()) && (rows == other.nrows()));

    Matrix<T> result(cols, rows, 0);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            if (data[i][k] == other[i][k])
                result[i][k] = 1;

    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::negate() &
{
    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            data[i][k] = -data[i][k];

    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::transpose() &
{
    Matrix<T> tmp(rows, cols);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            tmp[k][i] = data[i][k];

    swap(tmp);

    return *this;
}

template <typename T>
T Matrix<T>::trace() const
{
    assert(cols == rows);

    T result{};

    for (int i = 0; i < cols; ++i)
        result += data[i][i];

    return result;
}

template <typename T>
Matrix<T>::Matrix(int cols, int rows, T val) : data(new T *[rows]), cols(cols), rows(rows)
{
    for (int i = 0; i < rows; ++i)
    {
        data[i] = new T[cols];
        for (int k = 0; k < cols; ++k)
            data[i][k] = val;
    }
}

template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> initList) : Matrix<T>((*initList.begin()).size(), initList.size())
{
    int x_{}, y_;
    for (auto x : initList)
    {
        y_ = 0;
        for (auto y : x)
            data[x_][y_++] = y;
        ++x_;
    }
}

template <typename T>
void Matrix<T>::swap(Matrix<T> &rhs) noexcept
{
    std::swap(data, rhs.data);
    std::swap(cols, rhs.cols);
    std::swap(rows, rhs.rows);
}

template <typename T>
T **safe_copy(const T **src, int rows, int cols)
{
    T **dest = new T *[rows];
    for (int i = 0; i < rows; ++i)
        dest[i] = new T[cols];
    try
    {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                dest[i][j] = src[i][j];
    }
    catch (...)
    {
        for (int i = 0; i < rows; ++i)
        {
            delete[] dest[i];

            delete[] dest;
            throw;
        }
    }
    return dest;
}

template <typename T>
void Matrix<T>::dump(std::ostream &os) const
{
    for (int i = 0; i < rows; ++i)
    {
        for (int k = 0; k < cols; ++k)
            os << data[i][k] << ' ';
        os << std::endl;
    }
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &rhs) : data(safe_copy<T>((const T **)rhs.data, rows, cols)),
                                          rows(rhs.rows),
                                          cols(rhs.cols) {}

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs)
{
    Matrix tmp(rhs); // ex
                     //--------------------------------------//
    swap(tmp);       // noex
    return *this;
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> &&rhs) noexcept : data(rhs.data), rows(rhs.rows), cols(rhs.cols)
{
    rhs.data = nullptr;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&rhs) noexcept
{
    if (&rhs == this)
        return *this;

    swap(rhs);
    return *this;
}

template <typename T>
Matrix<T>::~Matrix()
{
    if (data != nullptr)
        for (int i = 0; i < rows; ++i)
            delete[] data[i];

    delete[] data;
}

#endif
