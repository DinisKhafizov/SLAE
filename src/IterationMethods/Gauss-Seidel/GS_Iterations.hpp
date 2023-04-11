#ifndef ITERS
#define ITERS

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> TopDownIteration(CSR A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x) {

    double Ux, Lx;
    int N = size(x), u;

    for (int k = 1; k < A.GetRow()[1]; ++k) {
		Ux += A.GetVal()[k] *  x[A.GetCol()[k]];
	}

    x[0] = (b[0] - Ux) / diag[0];

    for (int j = 1; j < N; ++j) {

        Ux = 0; //sum U * x_i
        Lx = 0; // sum L * x_(i+1)
        u = 1;

        for (int k = A.GetRow()[j]; k < A.GetRow()[j + 1] && A.GetCol()[k] < j; ++k) {
            Lx += A.GetVal()[k] * x[A.GetCol()[k]];
            u++; //counter of iterations (A - square matrix, so *num of counting Ux iterations*  = N - *number of counting Lx iterations*)
        }

        for (int k = A.GetRow()[j] + u; k < A.GetRow()[j + 1]; ++k) { //so here I use counter u
            Ux += A.GetVal()[k] * x[A.GetCol()[k]];
        }

        x[j] = (b[j] - Ux - Lx) / diag[j];

    }
    return x;
}

std::vector<double> DownUpIteration(CSR A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x) {
    double Ux, Lx;
    int N = size(x), u;
    for (int k = A.GetRow()[N - 1]; k < A.GetRow()[N] - 1; ++k) {
		Lx += A.GetVal()[k] * x[A.GetCol()[k]];
	}
    x.back() = (b[N - 1] - Lx) / diag[N - 1];

    for (int j = N - 1; j > 0; --j) {

        Ux = 0; //sum U * x_i
        Lx = 0; // sum L * x_(i+1)
        u = 1;

        for (int k = A.GetRow()[j - 1]; k < A.GetRow()[j] && A.GetCol()[k] < (j - 1); ++k) {
            Lx += A.GetVal()[k] * x[A.GetCol()[k]];
            u++; //counter of iterations (A - square matrix, so *num of counting Ux iterations*  = N - *number of counting Lx iterations*)
        }

        for (int k = A.GetRow()[j - 1] + u; k < A.GetRow()[j]; ++k) { //so here I use counter u
            Ux += A.GetVal()[k] * x[A.GetCol()[k]];
        }

        x[j - 1] = (b[j - 1] - Ux - Lx) / diag[j - 1];

    }
    return x;
}

#endif