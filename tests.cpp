#include <gtest/gtest.h>

#include "matrix.hpp"

TEST(test1, moveAssigment)
{
    Matrix matrix1(3, 2, 3);
    Matrix matrix2(3, 2, 1);

    matrix2 = std::move(matrix1);

    EXPECT_EQ(matrix2, Matrix(3, 2, 3));
}

TEST(test2, moveConstruct)
{
    Matrix matrix1(3, 2, 3);
    Matrix matrix2(std::move(matrix1));

    EXPECT_EQ(matrix2, Matrix(3, 2, 3));
}

TEST(test3, negate)
{
    Matrix matrix1(3, 2, 3);

    matrix1.negate();

    EXPECT_EQ(matrix1, Matrix(3, 2, -3));
}

TEST(test4, eye)
{
    Matrix<int> matrix1(3, 3);
    for (int i = 0; i < 3; ++i)
        matrix1[i][i] = 1;

    Matrix matrix2 = Matrix<int>::eye(3, 3);

    EXPECT_EQ(matrix1, matrix2);
}

TEST(test5, trace)
{
    Matrix matrix1(3, 3, 7);

    EXPECT_EQ(matrix1.trace(), 7 * 3);
}

TEST(test6, det)
{
    Matrix matrix1(4, 4, 0);

    int g = -1;

    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 4; ++k)
            matrix1[i][k] = g++;

    matrix1[2][0] = 10;
    matrix1[0][2] = 10;

    EXPECT_EQ(matrix1.det(), 432);
}

TEST(test7, equal)
{
    Matrix matrix1(4, 4, 5);

    EXPECT_EQ(matrix1.equal(matrix1), Matrix(4, 4, 1));
}

TEST(test8, less)
{
    Matrix matrix1(4, 4, 5);

    EXPECT_EQ(matrix1.less(matrix1), Matrix<int>(4, 4));
}

TEST(test9, transpose)
{
    Matrix<int> matrix1(4, 2);
    Matrix<int> matrix2(2, 4);

    matrix2.transpose();

    EXPECT_EQ(matrix1, matrix2);
}

TEST(test10, initList)
{
    Matrix matrix1(2, 2, 3);
    Matrix matrix2 = {{3, 3},
                      {3, 3}};

    EXPECT_EQ(matrix1, matrix2);
}

TEST(test11, types)
{
    Matrix matrixFloat = {{1.0f}};
    Matrix matrixBoolTrue = {{true}};
    Matrix matrixBoolFalse = {{false}};
    Matrix matrixMatrixInt(2, 2, Matrix<int>(2, 2, 1));

    EXPECT_TRUE     (matrixBoolTrue[0][0]);
    EXPECT_FALSE    (matrixBoolFalse[0][0]);
    EXPECT_FLOAT_EQ (matrixFloat[0][0], 1.0f);
    EXPECT_EQ       (matrixMatrixInt[0][0], Matrix<int>(2, 2, 1));
}


