#include "Acceleration/Cheb_Accel.hpp"
#include <gtest/gtest.h>

TEST(cheb, numb1) {
    std::vector<int> ns = {0, 1}, ns1;
    int n = 1;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ASSERT_NEAR(ns[i], ns1[i], 0.001);
    }
}
TEST(cheb, numb2) {
    std::vector<int> ns = {0, 3, 1, 2}, ns1;
    int n = 2;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ASSERT_NEAR(ns[i], ns1[i], 0.001);
    }
}
TEST(cheb, numb3) {
    std::vector<int> ns = {0, 7, 3, 4, 1, 6, 2, 5}, ns1;
    int n = 3;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ASSERT_NEAR(ns[i], ns1[i], 0.001);
    }
}
TEST(cheb, numb4) {
    std::vector<int> ns = {0, 15, 7, 8, 3, 12, 4, 11, 1, 14, 6, 9, 2, 13, 5, 10}, ns1;
    int n = 4;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ASSERT_NEAR(ns[i], ns1[i], 0.001);
    }
}
