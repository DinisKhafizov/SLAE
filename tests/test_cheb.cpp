#include "Acceleration/Cheb_Accel.hpp"
#include <gtest/gtest.h>

TEST(cheb, numb1) {
    std::vector<int> ns = {0, 1}, ns1;
    int n = 1;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ns[i] == ns1[i];
    }
}
TEST(cheb, numb2) {
    std::vector<int> ns = {0, 3, 1, 2}, ns1;
    int n = 2;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ns[i] == ns1[i];
    }
}
TEST(cheb, numb3) {
    std::vector<int> ns = {0, 7, 3, 4, 1, 6, 2, 5}, ns1;
    int n = 3;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ns[i] == ns1[i];
    }
}
TEST(cheb, numb4) {
    std::vector<int> ns = {0, 15, 7, 8, 3, 12, 4, 11, 1, 14, 6, 9, 2, 13, 5, 10}, ns1;
    int n = 4;
    ns1 = Numbers(n);
    for (int i = 0; i < size(ns); ++i) {
        ns[i] == ns1[i];
    }
}
TEST(cheb, taus1) {
    std::vector<double> taus = {0.03024603, 0.032245414, 0.036732, 0.0448933932, 0.05910861575, 
    0.08355057537, 0.122236866828, 0.16311062418441}, taus1;
    taus1 = TauFromCheb(8, 5.867, 33.326);
    for (int i = 0; i < size(taus1); ++i) {
        ASSERT_NEAR(taus1[i], taus[i], 0.000001);
    }
}