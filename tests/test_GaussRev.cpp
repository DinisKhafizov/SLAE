#include "Gauss/GaussReverse.hpp"
#include <gtest/gtest.h>



TEST(gauss, gauss_rev1) {
    std::vector<double> A = {1, 2, 7, 0, 4, 3, 0, 0, 10}, b = {3, 2, 1};
    std::vector<double> x, x1 = {1.45, 0.425, 0.1};
    Matrix B(A, 3, 3);
    x = GaussReverse(B, b);
    for(size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.01);
    }
}
TEST(gauss, gauss_rev2) {
    std::vector<double> A = {1.123, 2.123, 7.421, 3.521, 0, 45, 32.13, 241.2, 0, 0, 10, 312.55, 0, 0, 0, 32.1231}, b = {3, 2, 1, 7};
    std::vector<double> x, x1 = {39.4, 3.668, -6.711, 0.218};
    Matrix B(A, 4, 4);
    x = GaussReverse(B, b);
    for(size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.01);
    }  
}

TEST(gauss, gauss_rev3) {
    std::vector<double> A = {3}, b = {5};
    std::vector<double> x, x1 = {5./3};
    Matrix B(A, 1, 1);
    x = GaussReverse(B, b);
    for(size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.01);
    }  
}
TEST(gauss, gauss_rev4) {
    std::vector<double> A = {15, 13, 123, 1, 2, 0, 7, 113, 94, -24, 0, 0, 236, 85, 100, 0, 0, 0, -21, -123, 0, 0, 0, 0, 923};
    std::vector<double> b = {213, 23, 623, 531, 123};
    Matrix B(A, 5, 5);
    std::vector<double> x = GaussReverse(B, b), x1 = {-221.364, 160.519, 11.972, -26.066, 0.133};
    for(size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.01);
    }  
}
TEST(gauss, gauss_rev_spec1) {
    std::vector<double> A = {15, 13, 123, 1, 2, 0, 7, 113, 94, -24, 0, 0, 236, 85, 100, 0, 0, 0, -21, -123, 0, 0, 0, 0, 923};
    std::vector<double> b = {213, 23, 623, 531};
    Matrix B(A, 5, 5);
    std::vector<double> x = GaussReverse(B, b, 4), x1 = {-213.219, 153.207, 11.747, -25.286};
    for(size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.01);
    }  
}
TEST(gauss, gauss_rev_spec2) {
    std::vector<double> A = {1.123, 2.123, 7.421, 3.521, 0, 45, 32.13, 241.2, 0, 0, 10, 312.55, 0, 0, 0, 32.1231}, b = {3, 2, 1};
    std::vector<double> x, x1 = {2.062, -0.027, 0.1};
    Matrix B(A, 4, 4);
    x = GaussReverse(B, b, 3);
    for(size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.01);
    }  
}