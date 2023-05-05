#include "GS_Iterations.hpp"

std::vector<double> TopDownIteration(const CSR &A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x) {
    const int N = size(x);
    int row_begin, row_end, counter;
    for (size_t i = 0; i < N; ++i) {
        x[i] = b[i];
        counter = 1;
        row_begin = A.GetRow()[i];
        row_end = A.GetRow()[i + 1];
        for (size_t j = row_begin; j < row_end && A.GetCol()[j] < i; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
            ++counter;
        }
        for (size_t j = row_begin + counter; j < row_end; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
        } 
        x[i] /= diag[i];
    }
    return x;
}

std::vector<double> DownUpIteration(const CSR &A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x) {
    const int N = size(x);
    int row_begin, row_end, counter;

    for (int i = N - 1; i >= 0; --i) {
        x[i] = b[i];
        counter = 1;
        row_begin = A.GetRow()[i];
        row_end = A.GetRow()[i + 1];
        for (size_t j = row_begin; j < row_end && A.GetCol()[j] < i; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
            ++counter;
        }
        for (size_t j = row_begin + counter; j < row_end; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
        } 
        x[i] /= diag[i];
    }
    return x;
}