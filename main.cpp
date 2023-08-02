#include <iostream>
#include <memory>

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
    void show() const
    {
        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t k = 0; k < cols; ++k)
                std::cout << data[i][k] << ' ';
            std::cout << std::endl;
        }
    }

    bool equal(const Matrix &other) const;
    bool less(const Matrix &other) const;
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
};

template <typename T>
int Matrix<T>::ncols() const { return cols; }

template <typename T>
int Matrix<T>::nrows() const { return rows; }

template <typename T>
Matrix<T>::Matrix(int cols, int rows, T val) : data(new T *[rows]), cols(cols), rows(rows)
{
    for (size_t i = 0; i < rows; ++i)
    {
        data[i] = new T[cols];
        for (size_t k = 0; k < cols; ++k)
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
    for (size_t i = 0; i < rows; ++i)
        dest[i] = new T[cols];
    try
    {
        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < cols; ++j)
                dest[i][j] = src[i][j];
    }
    catch (...)
    {
        for (size_t i = 0; i < rows; ++i)
        {
            delete[] dest[i];

            delete[] dest;
            throw;
        }
    }
    return dest;
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
    for (size_t i = 0; i < rows; ++i)
        delete[] data[i];

    delete[] data;
}

int main(int argc, char const *argv[])
{
    Matrix<int> matrix1(3, 2, 3);
    Matrix matrix2(3, 2, 1);

    matrix2 = std::move(matrix1);

    matrix2.show();

    return 0;
}
