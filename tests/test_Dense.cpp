#include <gtest/gtest.h>
#include "Matrixes/SMatrix/StandardMatrix.hpp"

TEST(dense, transpose1) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 8, 2);
    Matrix B = A.transpose(1);
    for (size_t i = 0, end = A.GetM(); i < end; ++i) {
        for (size_t j = 0, end = A.GetN(); j < end; ++j) {
            ASSERT_DOUBLE_EQ(A(i, j), B(j, i));
        }
    }
}
TEST(dense, transpose2) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 4, 4);
    Matrix B = A.transpose(1);
    for (size_t i = 0, end = A.GetM(); i < end; ++i) {
        for (size_t j = 0, end = A.GetN(); j < end; ++j) {
            ASSERT_DOUBLE_EQ(A(i, j), B(j, i));
        }
    }
}
TEST(dense, transpose3) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 16, 1);
    Matrix B = A.transpose(1);
    for (size_t i = 0, end = A.GetM(); i < end; ++i) {
        for (size_t j = 0, end = A.GetN(); j < end; ++j) {
            ASSERT_DOUBLE_EQ(A(i, j), B(j, i));
        }
    }
}
TEST(dense, transpose4) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 1, 16);
    Matrix B = A.transpose(1);
    for (size_t i = 0, end = A.GetM(); i < end; ++i) {
        for (size_t j = 0, end = A.GetN(); j < end; ++j) {
            ASSERT_DOUBLE_EQ(A(i, j), B(j, i));
        }
    }
}
TEST(dense, mult1){
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 8, 2);
    std::vector<double> x = {5, 1}, res, res1 = {5, 38, 296.213, 114, 617, 25, 38, 123168};
    res = A * x;
    for(size_t i = 0, end = size(res); i < end; ++i) {
        ASSERT_NEAR(res[i], res1[i], 0.001);
    }
}
TEST(dense, partly_dot1) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    std::vector<double> x = {3, 4};
    Matrix A(a, 4, 4);
    std::vector<double> res = A.partly_dot(x),res1 = {3, 963.852, 377, 50};
    for(size_t i = 0, end = size(res); i < end; ++i) {
        ASSERT_NEAR(res[i], res1[i], 0.001);
    }
}

TEST(dense, dot_matrix1) {
    std::vector<double> a = {1, 2, 2, 3, 4, 5, 7, 8, 11}, b = {7, 6, 1, 3, 5 ,8, 8 , 9, 13, 23, 32, 10}, c = {43, 68, 81, 41, 106, 165, 195, 95, 232, 359, 423, 203};
    Matrix A(a, 3, 3), B(b, 3, 4), res, res1(c, 3, 4);
    res = A * B;
    for (size_t i = 0, end = A.GetM(); i < end; ++i) {
        for (size_t j = 0, end = A.GetN(); j < end; ++j) {
            ASSERT_DOUBLE_EQ(res(i, j), res1(i, j));
        }
    }
}
