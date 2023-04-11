#ifndef GAUSS_SEID_METH
#define GAUSS_SEID_METH

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> GSM(CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {

	int N = A.GetN(), u; 

	double Ux = 0, Lx = 0; 
	
	std::vector<double> x = x_0, vec(N), x_last(N, 100), diag(N);

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
			}
		}
	}

	while (first_norm(x - x_last) > tolerance) {

		x_last = x;

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

	}

    return x;
}


#endif