#include <gtest/gtest.h>
#include "Matrixes/SMatrix/StandardMatrix.hpp"

TEST(dense, transpose1) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 8, 2);
    A.show();
    A.transpose();
    A.show();
}
TEST(dense, transpose2) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 4, 4);
    A.show();
    A.transpose();
    A.show();
}
TEST(dense, transpose3) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 16, 1);
    A.show();
    A.transpose();
    A.show();
}
TEST(dense, transpose4) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 1, 16);
    A.show();
    A.transpose();
    A.show();
}
TEST(dense, mult1){
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    Matrix A(a, 8, 2);
    std::vector<double> x = {5, 1}, res;
    res = A * x;
}
TEST(dense, partly_dot1) {
    std::vector<double> a = {1, 0, 7, 3, 13, 231.213, 21, 9, 123, 2, 4, 5, 6, 8, 9, 123123};
    std::vector<double> x = {3, 4};
    Matrix A(a, 4, 4);
    A.show();
    std::vector<double> res = A.partly_dot(x);
    for (size_t i = 0; i < 4; ++i) {
        std::cout << res[i] << std::endl;
    }
}