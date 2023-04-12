#include "Vect/VectorOperations.hpp"
#include <gtest/gtest.h>
#include"Matrixes/CSR/MatrixOnCSR.hpp"

TEST(mult, mult1) {
    std::vector<double> a = {3, 5, 1}, b = {2, 7, 1};
    double res = 314, res1;
    std::vector<double> vals = {1, 2, 3, 4, 5, 3, 1, 2, 1};
    std::vector<int> cols = {0, 1, 2, 0, 1, 2, 0, 1, 2}, rows = {0, 3, 6, 9};
    CSR A(vals, cols, rows, 3);
    double x = 2;
    res1 = a * (A * b);
    ASSERT_DOUBLE_EQ(res, res1);
}

/*
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/