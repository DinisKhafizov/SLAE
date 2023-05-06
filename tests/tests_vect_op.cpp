#include "Vect/VectorOperations.hpp"
#include <gtest/gtest.h>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Matrixes/SMatrix/StandardMatrix.hpp"


TEST(vec, mult2) {
    std::vector<double> x = {3, 2, 1}, vals = {3, 7, 5, 4, 8, 9, 13, 15, 17};
    std::vector<int> cols = {0, 1, 2, 0, 1, 2, 0, 1, 2}, rows = {0, 3, 6, 9};
    CSR A(vals, cols, rows, 3);
    x = A*x;
    for (size_t i = 0; i < 3; ++i){
        std::cout << x[i] << std::endl;
    }
}
TEST(vec, dot1) {
    std::vector<double> a = {7, 8, 1, 5}, b = {2, 3, 1};
    Matrix res = dot(a, b);
    res.show();
}
TEST(vec, self_dot1) {
    std::vector<double> a = {7, 8, 1, 5}, b = {2, 3, 1};
    Matrix res = dot_self(a);
    res.show();
}

/*
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/