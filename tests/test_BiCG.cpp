#include "ProjectionMethods/BiCG/BiCG.hpp"
#include <gtest/gtest.h>

TEST(Bicg, solve1) {
    std::vector<double> b = {12, 51, 3, 123}, x_0={0, 0, 0, 0}, vals = {5, -1, -2, -1, 4, 1, 1, 3, -1, -2, -1, 4}, res = {23.718, 15.987, 10.772, 45.302}, res1;
    std::vector<int> rows = {0, 3, 6, 9, 12}, cols = {0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3};
    CSR A(vals, cols, rows, 4);
    res1 = BiCG(A, b, x_0, 0.000001);
    for (int i = 0; i < 4; ++i) {
        ASSERT_NEAR(res[i], res1[i], 0.001);
    }
}