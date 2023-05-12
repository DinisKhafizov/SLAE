#include <gtest/gtest.h>
#include"Matrixes/CSR/MatrixOnCSR.hpp"


TEST(csr, op_skobki1) {
    std::vector<double> values = {1.};
    std::vector<int> cols = {0};
    std::vector<int> rows = {0, 1};
    std::vector<double> check = {1.};
    int Ncols = 1;
    CSR A(values, cols, rows, Ncols);
    for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < 1; ++j) {
           ASSERT_DOUBLE_EQ(check[i*3 + j], A(i, j));
        }
    }
}


TEST(csr, op_skobki3) {
    std::vector<double> values = {5., 3., 2., 7., 4., 6.};
    std::vector<double> check = {5., 3., 2., 7., 4., 6.};
    std::vector<int> cols = {0, 1, 2, 3, 4, 5};
    std::vector<int> rows = {0, 6};
    int Ncols = 6;
    CSR A(values, cols, rows, Ncols);
    for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < 6; ++j) {
           ASSERT_DOUBLE_EQ(check[i*3 + j], A(i, j));
        }
    }
}
TEST(csr, op_skobki4) {
    std::vector<double> values = {5., 3., 2., 7., 4., 6.};
    std::vector<double> check = {5., 3., 2., 7., 4., 6.};
    std::vector<int> cols = {0, 0, 0, 0, 0, 0};
    std::vector<int> rows = {0, 1, 2, 3, 4, 5, 6};
    int Ncols = 3;
    CSR A(values, cols, rows, Ncols);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 1; ++j) {
           ASSERT_DOUBLE_EQ(check[i], A(i, j));
        }
    }
}

TEST(csr, op_skobki5) {
    std::vector<double> values = {5., 3., 2., 7., 4., 6.};
    std::vector<int> cols = {0, 1, 0, 2, 0, 1};
    std::vector<int> rows = {0, 2, 4, 6};
    std::vector<double> check = {5, 3, 0, 2, 0, 7, 4, 6, 0};
    int Ncols = 3;
    CSR A(values, cols, rows, Ncols);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
           ASSERT_DOUBLE_EQ(check[i*3 + j], A(i, j));
        }
    }
}

TEST(csr, matrix_multi_vec1) {
    std::vector<double> values = {5., 3., 2., 7., 4., 6.};
    std::vector<int> cols = {0, 1, 0, 2, 0, 1};
    std::vector<int> rows = {0, 2, 4, 6};
    std::vector<double> f = {3., 5., 9.};
    std::vector<double> check = {30., 69., 42.}, res;
    int Ncols = 3;
    CSR A(values, cols, rows, Ncols);
    res = A * f;
    for (int i = 0; i < 3; ++i) {
        ASSERT_DOUBLE_EQ(res[i], check[i]);
    }
}

TEST(csr, matrix_multi_vec2) {
    std::vector<double> values = {10., 200., 35., 43., 52., 61., 47., 38., 59., 110., 511., 1200., 135., 141., 15.3};
    std::vector<int> cols = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    std::vector<int> rows = {0, 5, 10, 15};
    std::vector<double> f = {72., 54., 99., 12.123, 3};
    std::vector<double> check = {15662.289, 11737.257, 116712.243}, res;
    int Ncols = 5;
    CSR A(values, cols, rows, Ncols);
    res = A * f;
    for (int i = 0; i < 3; ++i) {
        ASSERT_DOUBLE_EQ(res[i], check[i]);
    }
}

TEST(csr, matrix_multi_vec3) {
    std::vector<double> values = {5., 3., 2., 7., 4.};
    std::vector<int> cols = {0, 1, 2, 3, 4};
    std::vector<int> rows = {0, 5};
    std::vector<double> f = {23., 35.21, 123.99, 213.23, 864.23};
    std::vector<double> check = {5418.14}, res;
    int Ncols = 3;
    CSR A(values, cols, rows, Ncols);
    res = A * f;
    for (int i = 0; i < 1; ++i) {
        ASSERT_DOUBLE_EQ(res[i], check[i]);
    }
}

TEST(csr, csr_show1) {
    std::vector<double> values = {5., 3., 2., 7., 4.};
    std::vector<int> cols = {0, 1, 2, 3, 4};
    std::vector<int> rows = {0, 5};
    std::vector<double> f = {23., 35.21, 123.99, 213.23, 864.23};
    std::vector<double> check = {5418.14}, res;
    int Ncols = 5;
    CSR A(values, cols, rows, Ncols);
}
TEST(csr, csr_show2) {
    std::vector<double> values = {10., 200., 35., 43., 52., 61., 47., 38., 59., 110., 511., 1200., 135., 141., 15.3};
    std::vector<int> cols = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    std::vector<int> rows = {0, 5, 10, 15};
    std::vector<double> f = {72., 54., 99., 12.123, 3};
    std::vector<double> check = {15662.289, 11737.257, 116712.243}, res;
    CSR A(values, cols, rows, 7);
}
TEST(csr, transp1) {
    std::vector<double> values = {10., 200., 35., 43., 52., 61., 47., 38., 59., 110., 511., 1200., 135., 141., 15.3};
    std::vector<int> cols = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    std::vector<int> rows = {0, 5, 10, 15};
    std::vector<double> f = {72., 54., 99., 12.123, 3};
    std::vector<double> check = {15662.289, 11737.257, 116712.243}, res;
    CSR A(values, cols, rows, 7);
    A.transpose();
}