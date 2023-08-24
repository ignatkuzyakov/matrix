#include <iostream>
#include <memory>
#include <gtest/gtest.h>


template <typename T>
class Matrix
{
    int rows, cols;
    T **data = nullptr;

public:
    Matrix(const Matrix &rhs);
    Matrix(Matrix &&rhs);

    Matrix &operator=(const Matrix &rhs);
    Matrix &operator=(Matrix &&rhs);

    ~Matrix();

public:
    Matrix(int cols, int rows, T val = T{});

    template <typename It>
    Matrix(int cols, int rows, It start, It fin);

    static Matrix eye(int n, int m);

public: // operations
    Matrix &negate() &;
    Matrix &transpose() &;

public:
    int ncols() const;
    int nrows() const;

    T trace() const;

    void swap(Matrix &rhs) noexcept;

    bool operator==(const Matrix<T> &rhs) const ;

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
    Matrix<T> result(n, m, 0);
    if (n == m)

        for (int i = 0; i < ncols; ++i)
            result[i][i] = 1;

    else
        std::cout << "cols != rows" << std::endl;

    return result;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T> &rhs) const
{
    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            if (data[i][k] != rhs[i][k])
                return false;
                
    return true;
}

template <typename T>
Matrix<T> Matrix<T>::less(const Matrix &other) const
{
    Matrix<T> result(ncols, nrows, 0);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            if (data[i][k] < other[i][k])
                result[i][k] = 1;

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::equal(const Matrix &other) const
{
    Matrix<T> result(ncols, nrows, 0);

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
    Matrix<T> tmp(nrows, ncols);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
            tmp[k][i] = data[i][k];

    swap(tmp);

    return *this;
}

template <typename T>
T Matrix<T>::trace() const
{
    T result{};

    if (ncols == nrows)
        for (int i = 0; i < ncols; ++i)
            result += data[i][i];
    else
        std::cout << "cols != rows" << std::endl;

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
    if (rhs.cols != cols || rhs.rows != rows)
    {
        std::cout << "rhs.cols != cols || rhs.rows != rows" << std::endl;
        return *this;
    }
    Matrix tmp(rhs); // ex
                     //--------------------------------------//
    swap(tmp);       // noex
    return *this;
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> &&rhs) : data(rhs.data), rows(rhs.rows), cols(rhs.cols)
{
    rhs.data = nullptr;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&rhs)
{
    if (&rhs == this)
        return *this;

    swap(rhs);
    return *this;
}

template <typename T>
Matrix<T>::~Matrix()
{
    for (int i = 0; i < rows; ++i)
        delete[] data[i];

    delete[] data;
}


TEST(test1, matrix_move)
{
    Matrix<int> matrix1(3, 2, 3);
    Matrix<int> matrix2(3, 2, 1);
    Matrix<int> matrix3(3, 2, 3);

    matrix2 = std::move(matrix1);

    EXPECT_EQ(matrix2[0][0], matrix3[0][0]);
}

TEST(test2, negate)
{
    Matrix<int> matrix1(3, 2, 3);
    Matrix<int> matrix2(3, 2, -3);

    matrix1.negate();

    EXPECT_EQ(matrix1, matrix2);  
}



