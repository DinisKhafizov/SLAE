#include "Gauss/GaussReverse.hpp"
#include <gtest/gtest.h>

TEST(gauss, gauss_rev1) {
    std::vector<double> A = {1, 2, 7, 0, 4, 3, 0, 0, 10}, b = {3, 2, 1};
    std::vector<double> x;
    Matrix B(A, 3, 3);
    x = GaussReverse(B, b);
    for (int i = 0; i < 3; ++i) {
        std::cout << x[i] << std::endl;
    }
}
TEST(gauss, gauss_rev2) {
    std::vector<double> A = {1.123, 2.123, 7.421, 3.521, 0, 45, 32.13, 241.2, 0, 0, 10, 312.55, 0, 0, 0, 32.1231}, b = {3, 2, 1, 7};
    std::vector<double> x;
    Matrix B(A, 4, 4);
    x = GaussReverse(B, b);
    for (int i = 0; i < 4; ++i) {
        std::cout << x[i] << std::endl;
    }    
}