#include "Vect/VectorOperations.hpp"
#include <gtest/gtest.h>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Matrixes/SMatrix/StandardMatrix.hpp"


TEST(vec, mult2) {
    std::vector<double> x = {3, 2, 1}, vals = {3, 7, 5, 4, 8, 9, 13, 15, 17}, x1 = {28, 37, 86};
    std::vector<int> cols = {0, 1, 2, 0, 1, 2, 0, 1, 2}, rows = {0, 3, 6, 9};
    CSR A(vals, cols, rows, 3);
    x = A*x;
    for (size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_DOUBLE_EQ(x[i], x1[i]);
    }
}
TEST(vec, dot1) {
    std::vector<double> a = {7, 8, 1, 5}, b = {2, 3, 1};
    std::vector<double> vals = {14, 21, 7, 16, 24, 8, 2, 3, 1, 10, 15, 5};
    Matrix res1(vals, 4, 3);
    Matrix res = dot(a, b);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 3; ++j)
            ASSERT_DOUBLE_EQ(res(i, j), res1(i, j));
    }
}
TEST(vec, self_dot1) {
    std::vector<double> a = {7, 8, 1, 5};
    Matrix res = dot_self(a);
    std::vector<double> vals = {49, 56, 7, 35, 56, 64, 8, 40, 7, 8, 1, 5, 35, 40, 5, 25};
    Matrix res1(vals, 4, 4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j)
            ASSERT_DOUBLE_EQ(res(i, j), res1(i, j));
    }
}
