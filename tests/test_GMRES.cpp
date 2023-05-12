#include <gtest/gtest.h>
#include "ProjectionMethods/GMRES/GMRES.hpp"

//need to test: constr, newiter, givens(last_iter), newiter + givens(last_iter), getting V without cols

TEST(hessenberg, constructor1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Hessenberg heis(A, r_0);
}
TEST(hessenberg, constructor2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0);
}
TEST(hessenberg, constructor3) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0, 0);
}
TEST(hessenberg, newiter1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Hessenberg heis(A, r_0);
    for (size_t i = 0; i < 3; ++i) {
        heis.newIter(A);
    }
}
TEST(hessenberg, newiter2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0);
    for (size_t i = 0; i < 2; ++i) {
        heis.newIter(A);
    }
}
TEST(hessenberg, newiter3) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0, 0);
    for (size_t i = 0; i < 3; ++i) {
        heis.newIter(A);
    }
}
TEST(hessenberg, givens1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Hessenberg heis(A, r_0, 4);
    heis.givens();
}
TEST(hessenberg, iter_plus_givens1) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11}, r_0 = {1, 1, 1, 1};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9};
    CSR A(vals, cols, rows, 4);
    Hessenberg heis(A, r_0);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
} 
TEST(hessenberg, givens2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0, 3);
    heis.givens();
    
}
TEST(hessenberg, iter_plus_givens2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
    heis.newIter(A);
    heis.givens_last_iter();
}
TEST(hessenberg, V_getting1) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0);
}
TEST(hessenberg, V_getting2) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0, 2);
    Matrix V_ex_last = heis.get_V_exc_lastcol();
}
TEST(hessenberg, V_getting3) {
    std::vector<double> vals = {1, 2, 3, 4, 5, 6, 7, 8, 10}, r_0 = {1, 1, 1};
    std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> b = {1, 2, 9};
    CSR A(vals, cols, rows, 3);
    Hessenberg heis(A, r_0, 3);
    Matrix V_ex_last = heis.get_V_exc_lastcol();
}
TEST(gmres, slae1_recursive) {
    std::vector<double> vals = {4, 1, 3, 2, 1, 5, 6, 7, 3, 6, 9, 8, 2, 7, 8, 11};
    std::vector<int> rows = {0, 4, 8, 12, 16}, cols = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    std::vector<double> b = {1, 2, 3, 9}, x, x_0 = {0, 0, 0, 0}, x1 = {-5.182, -22.455, 7.818, 10.364};
    CSR A(vals, cols, rows, 4);
    x = GMRES(A, b, x_0, 0.00001);
    for (size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.001);
    }
}
TEST(gmres, slae2_recursive) {
    std::vector<double> vals = {3, 1, 2, 1, 1, 4, 1, 3, 2, 1, 6, 1, 2, 3, 1, 5, 2, 1, 2, 2, 4};
    std::vector<int> rows = {0, 4, 9, 12, 17, 21}, cols = {0, 1, 3, 4, 0, 1, 2, 3, 4, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 3, 4};
    std::vector<double> b = {1, 2, 3, 9, 0.5}, x, x_0 = {0, 0, 0, 0, 0}, x1 = {-1.24, -1.516, 0.192, 3.362, -0.488};
    CSR A(vals, cols, rows, 5);
    x = GMRES(A, b, x_0, 0.00001);
    for (size_t i = 0, end = size(x); i < end; ++i) {
        ASSERT_NEAR(x[i], x1[i], 0.001);
    }
}
