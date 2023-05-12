#include <gtest/gtest.h>
#include "LU/LU/ILU/ILU.hpp"

TEST(lu, ichol1) {
    std::vector<double> vals = {6, 1, -1, 1, 7, -1, -1, 8, 1, -1, -1, 1, 9, 2, -1, 2, 10};
    std::vector<int> cols = {0, 1, 2, 0, 1, 3, 0, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4}, rows = {0, 3, 6, 10, 14, 17};
    CSR A(vals, cols, rows, 5);
    Matrix L = Ichol_null(A), L_t;
    L_t = L;
}