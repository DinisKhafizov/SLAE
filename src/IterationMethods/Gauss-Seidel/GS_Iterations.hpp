#ifndef ITERS
#define ITERS

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

//Комменты:

//11.04: улучшил итерации сверху вниз и наоборот. Лишнюю итерацию с прибавлением x[i] += x[i] * diag[i] убрал
//посредством зануления диагональных элементов у матрицы А (у вектора values CSR матрицы). Других идей по оптимизации итераций пока нет.

std::vector<double> TopDownIteration(CSR A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x) {
    const int N = size(x);
    int r;
    for (int i = 0; i < N; ++i) {
        x[i] = b[i];
        r = A.GetRow()[i + 1];
        for (int j = A.GetRow()[i]; j < r; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
        }
        x[i] /= diag[i];
    }
    return x;
}

std::vector<double> DownUpIteration(CSR A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x) {
    const int N = size(x);
    int r;
    for (int i = N - 1; i >= 0; --i) {
        x[i] = b[i];
        r = A.GetRow()[i + 1];
        for (int j = A.GetRow()[i]; j < r; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
        }
        x[i] /= diag[i];
    }
    return x;
}

#endif