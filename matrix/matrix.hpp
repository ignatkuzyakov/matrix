#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

template <typename T>
class MatrixBuf
{
protected:
    int rows, cols;
    T **data;

protected:
    MatrixBuf &operator=(const MatrixBuf &rhs) = delete;
    MatrixBuf &operator=(MatrixBuf &&rhs) noexcept;

    MatrixBuf(const MatrixBuf &rhs) = delete;
    MatrixBuf(MatrixBuf &&rhs) noexcept;

    MatrixBuf(int cols, int rows);
    MatrixBuf(T **data, int cols, int rows);
    ~MatrixBuf();

    void swap(MatrixBuf &rhs) noexcept;
};

template <typename T>
MatrixBuf<T>::MatrixBuf(int cols, int rows) : data((rows || cols) ? (new T *[rows]) : (nullptr)), cols(cols), rows(rows)
{
    int constructed = 0;
    try
    {
        for (; constructed < rows; ++constructed)
            data[constructed] = new T[cols];
    }
    catch (...)
    {
        for (int k = 0; k < constructed; ++k)
            delete data[k];
        throw;
    }
}

template <typename T>
MatrixBuf<T>::MatrixBuf(T **data, int cols, int rows) : data(data), cols(cols), rows(rows) {}

template <typename T>
MatrixBuf<T>::MatrixBuf(MatrixBuf<T> &&rhs) noexcept : data(rhs.data), rows(rhs.rows), cols(rhs.cols)
{
    rhs.data = nullptr;
    rhs.rows = 0;
    rhs.cols = 0;
}

template <typename T>
MatrixBuf<T> &MatrixBuf<T>::operator=(MatrixBuf<T> &&rhs) noexcept
{
    swap(rhs);
    return *this;
}

template <typename T>
MatrixBuf<T>::~MatrixBuf()
{
    for (int i = 0; i < rows; ++i)
        delete[] data[i];

    delete[] data;
}

template <typename T>
void MatrixBuf<T>::swap(MatrixBuf<T> &rhs) noexcept
{
    std::swap(data, rhs.data);
    std::swap(cols, rhs.cols);
    std::swap(rows, rhs.rows);
}

template <typename T>
class Matrix : private MatrixBuf<T>
{
    using MatrixBuf<T>::data;
    using MatrixBuf<T>::rows;
    using MatrixBuf<T>::cols;
    using MatrixBuf<T>::swap;

public:
    Matrix &operator=(const Matrix &rhs);
    Matrix &operator=(Matrix &&rhs) noexcept = default;

    Matrix(const Matrix &rhs);
    Matrix(Matrix &&rhs) noexcept = default;

    Matrix(std::initializer_list<std::initializer_list<T>> list);
    explicit Matrix(int cols = 0, int rows = 0, T val = {});

public:
    static Matrix eye(int n, int m);

public: // operations
    Matrix &negate() &;
    Matrix &transpose() &;
    Matrix &invert() &;

public:
    int ncols() const;
    int nrows() const;

    T trace() const;

    // The det method using the sum of the products
    // of the elements of any one row or column and their cofactors.
    T det() const;

private:
    T det(Matrix matrix) const;

public:
    bool operator==(const Matrix &rhs) const;

    Matrix operator+(const Matrix &rhs) const;
    Matrix &operator+=(const Matrix &rhs);
    Matrix operator+(const T &val) const;
    Matrix &operator+=(const T &val);
    Matrix operator++(int);
    Matrix &operator++();
    Matrix operator+() const;

    Matrix operator-(const Matrix &rhs) const;
    Matrix &operator-=(const Matrix &rhs);
    Matrix operator-(const T &val) const;
    Matrix &operator-=(const T &val);
    Matrix operator--(int);
    Matrix &operator--();
    Matrix operator-() const;

    Matrix operator*(const Matrix &rhs) const;
    Matrix &operator*=(const Matrix &rhs);
    Matrix operator*(const T &val) const;
    Matrix &operator*=(const T &val);

    Matrix operator/(const Matrix &rhs) const;
    Matrix &operator/=(const Matrix &rhs);
    Matrix operator/(const T &val) const;
    Matrix &operator/=(const T &val);

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
T **safe_copy(T **src, int rows, int cols)
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
Matrix<T>::Matrix(const Matrix<T> &rhs) : MatrixBuf<T>(safe_copy<T>(rhs.data, rhs.rows, rhs.cols), rhs.cols, rhs.rows) {}

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs)
{
    Matrix tmp(rhs); // ex
                     //--------------------------------------//
    swap(tmp);       // noex
    return *this;
}

template <typename T>
Matrix<T>::Matrix(int cols, int rows, T val) : MatrixBuf<T>(cols, rows)
{
   Matrix dest (*this);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            dest[i][k] = val;

    swap(dest);
}

template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> initList) : MatrixBuf<T>((*initList.begin()).size(), initList.size())
{
    Matrix<T> dest (*this);

    int rws = 0, cls;
    for (auto x : initList)
    {
        cls = 0;
        for (auto y : x)
            dest[rws][cls++] = y;
        ++rws;
    }

    swap(dest);
}

template <typename T>
Matrix<T> &Matrix<T>::invert() &
{
    T vl = this->det();

    assert(vl != 0);

    Matrix<T> subMatrix(cols - 1, cols - 1);

    Matrix tmp(*this);

    for (int i = 0; i < cols; i++)
    {
        for (int k = 0; k < cols; k++)
        {
            for (int g = 0; g < cols; g++)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (j == k || i == g)        continue;
                    else if ((j < k) && (g > i)) subMatrix[g - 1][j] = data[g][j];
                    else if ((j > k) && (g > i)) subMatrix[g - 1][j - 1] = data[g][j];
                    else if ((j > k) && (g < i)) subMatrix[g][j - 1] = data[g][j];
                    else                         subMatrix[g][j] = data[g][j];
                }
            }
            tmp[i][k] = std::pow(-1, k + i) * subMatrix.det();
        }
    }

    tmp.transpose();
    tmp  /= vl;
    swap(tmp);

    return *this;
}

template <typename T>
T Matrix<T>::det() const
{
    assert(cols == rows);
    return det(*this);
}

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
                    if (j == k)     continue;
                    else if (j < k) subMatrix[i - 1][j] = matrix[i][j];
                    else            subMatrix[i - 1][j - 1] = matrix[i][j];
                }
            }
            resultD += std::pow(-1, k + 2) * matrix[0][k] * det(subMatrix);
        }
        return resultD;
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const
{
    assert((rhs.ncols() == cols) && (rhs.nrows() == rows));

    Matrix<T> tmp(rhs);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            tmp[i][j] += data[i][j];

    return tmp;
}

template <typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &rhs) { return (*this = rhs + *this); }

template <typename T>
Matrix<T> Matrix<T>::operator+() const { return *this; }

template <typename T>
Matrix<T> Matrix<T>::operator-() const
{
    Matrix<T> tmp(*this);

    return tmp.negate();
}

template <typename T>
Matrix<T> Matrix<T>::operator++(int)
{
    Matrix<T> tmp(*this);

    *this += Matrix<T>(cols, rows, 1);

    return tmp;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const T &val) const { return (*this + Matrix<T>(cols, rows, val)); }

template <typename T>
Matrix<T> &Matrix<T>::operator+=(const T &val) { return (*this = *this + val); }

template <typename T>
Matrix<T> &Matrix<T>::operator++() { return (*this += Matrix<T>(cols, rows, 1)); }

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const
{
    Matrix<T> tmp(*this);

    return (tmp += (-rhs));
}

template <typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &rhs) { return (*this = *this - rhs); }

template <typename T>
Matrix<T> Matrix<T>::operator-(const T &val) const { return (*this - Matrix<T>(cols, rows, val)); }

template <typename T>
Matrix<T> &Matrix<T>::operator-=(const T &val) { return (*this = *this - val); }

template <typename T>
Matrix<T> Matrix<T>::operator--(int)
{
    Matrix<T> tmp(*this);

    *this -= Matrix<T>(cols, rows, 1);

    return tmp;
}

template <typename T>
Matrix<T> &Matrix<T>::operator--() { return (*this -= Matrix<T>(cols, rows, 1)); }

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const
{
    assert(cols == rhs.nrows());

    Matrix<T> tmp(rhs.cols, rows);

    for (int i = 0, nrows = tmp.nrows(); i < nrows; ++i)
        for (int j = 0, ncols = tmp.ncols(); j < ncols; ++j)
            for (int g = 0; g < cols; ++g)
                tmp[i][j] += (data[i][g] * rhs[g][j]);

    return tmp;
}

template <typename T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &rhs) { return (*this = *this * rhs); }

template <typename T>
Matrix<T> Matrix<T>::operator*(const T &val) const
{
    Matrix<T> tmp(*this);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            tmp[i][j] *= val;

    return tmp;
}

template <typename T>
Matrix<T> &Matrix<T>::operator*=(const T &val) { return (*this = *this * val); }

template <typename T>
Matrix<T> Matrix<T>::operator/(const Matrix<T> &rhs) const
{
    Matrix tmp(rhs);

    return (*this * tmp.invert());
}

template <typename T>
Matrix<T> &Matrix<T>::operator/=(const Matrix<T> &rhs) { return (*this = (*this / rhs)); }

template <typename T>
Matrix<T> Matrix<T>::operator/(const T &val) const
{
    Matrix<T> tmp(*this);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            tmp[i][j] /= val;

    return tmp;
}

template <typename T>
Matrix<T> &Matrix<T>::operator/=(const T &val) { return (*this = *this / val); }

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
void Matrix<T>::dump(std::ostream &os) const
{
    for (int i = 0; i < rows; ++i)
    {
        for (int k = 0; k < cols; ++k)
            os << data[i][k] << ' ';
        os << std::endl;
    }
}

#endif
