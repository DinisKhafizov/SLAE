#include "ILU.hpp"

/*
std::pair<Matrix, Matrix> ILU_dec(Matrix A, const int p) {
    const int N = A.GetN();
    Matrix L(N, N), lvl(N, N);
    double new_lvl;

    for (size_t i = 0; i< N; ++i) {
        L(i, i) = 1;
        for (size_t j = 0; j < N; ++j) {
            if (A(i, j) != 0) {
                lvl(i, j) = 0;
            }
            else {
                lvl(i, j) = N - 1;
            }
        }
    }
    for (size_t i = 0; i < N; ++i) {
        for (int k = 0, end = i - 1; k < end; ++k) {
            if ((lvl(i, k) <= p) && (A(k, k) != 0)) {
                A(i, k) = A(i, k) - A(k, k); //sign
            }
            for (size_t j = k; j < N; ++j) {
                A(i, j) = A(i, j) - A(i, k) * A(k, j);
                new_lvl = lvl(i, k) + lvl(k, j) + 1;
                if (lvl(i,j) > new_lvl) {
                    lvl(i, j) = new_lvl;
                }
            }
        }
        for (size_t k = 0; k < N; ++k) {
            if (lvl(i, k) > p) {
                A(i, k) = 0;
            }
        }
    }
    std::pair<Matrix, Matrix> LU;
    LU.first = L;
    LU.second = A;
    lvl.show();
    return LU;
}

*/

Matrix Ichol_null(const CSR &A) {
    const int N = A.GetN();
    double sum, a_ij;
    Matrix L(N, N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            a_ij = A(i, j);
            if (a_ij != 0) {
                sum = L.getRow(i, j) * L.getRow(j, j);
                if (i == j) {
                    L(i, j) = sqrt(std::abs(a_ij - sum));
                }
                else{
                    L(i, j) = (a_ij - sum)/L(j, j);
                }
            }
        }
    }
    return L;
}