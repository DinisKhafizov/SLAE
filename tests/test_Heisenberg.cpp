#include <gtest/gtest.h>
#include "ProjectionMethods/ArnoldiAlgorythm/ArnoldiAlg.hpp"

//need to test: constr, newiter, givens(last_iter), newiter + givens(last_iter), getting V without cols

TEST(heisenberg, constructor1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Heisenberg heis(A, r_0);
    heis.get_H().show();
    heis.get_V().show();
}
TEST(heisenberg, constructor2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0);
    heis.get_H().show();
    heis.get_V().show();
}
TEST(heisenberg, constructor3) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0, 0);
    heis.get_H().show();
    heis.get_V().show();
}
TEST(heisenberg, newiter1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Heisenberg heis(A, r_0);
    for (size_t i = 0; i < 3; ++i) {
        heis.newIter(A);
        heis.get_H().show();
        heis.get_V().show();
    }
}
TEST(heisenberg, newiter2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0);
    for (size_t i = 0; i < 2; ++i) {
        heis.newIter(A);
        heis.get_H().show();
        heis.get_V().show();
    }
}
TEST(heisenberg, givens1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Heisenberg heis(A, r_0, 4);
    heis.givens();
    heis.get_H().show();
    heis.get_V().show();
}
TEST(heisenberg, iter_plus_givens1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Heisenberg heis(A, r_0);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.get_H().show();
    heis.get_V().show();
} 
TEST(heisenberg, givens2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0, 3);
    heis.givens();
    heis.get_H().show();
    heis.get_V().show();
    
}
TEST(heisenberg, iter_plus_givens2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.get_H().show();
    heis.get_V().show();
}
TEST(heisenberg, V_getting1) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0);
    heis.get_V_exc_lastcol().show();
    heis.get_V().show();
}
TEST(heisenberg, V_getting2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0, 2);
    heis.get_V_exc_lastcol().show();
    heis.get_V().show();
}
TEST(heisenberg, V_getting3) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Heisenberg heis(A, r_0, 3);
    heis.get_V_exc_lastcol().show();
    heis.get_V().show();
}