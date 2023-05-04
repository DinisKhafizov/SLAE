#include <gtest/gtest.h>
#include "ProjectionMethods/ArnoldiAlgorythm/ArnoldiAlg.hpp"
#include "Gauss/GaussReverse.hpp"

TEST(heisenberg, arnoldi_1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 3);
    std::vector<double> res, x_0 = {0, 0, 0, 0};
    res = GMRES(A, b, x_0, 0.000001, 3);

}