#include <gtest/gtest.h>
#include "ProjectionMethods/ArnoldiAlgorythm/ArnoldiAlg.hpp"

TEST(heisenberg, arnoldi_1) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    double k;
    std::vector<double> kek;
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0, 5);
    Matrix H = heis.get_H(), V = heis.get_V();
    H.show();
    V.show();

}