#ifndef GAUSS_SEID_METH
#define GAUSS_SEID_METH

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> GSM(CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {
	const int N = size(x_0);
	std::vector<double> x = x_0, x_last(N, 100), diag(N);

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
				A.GetVal()[k] = 0;
			}
		}
	}
	while (first_norm(x - x_last) > tolerance) {
		x_last = x;
    	for (int i = 0; i < N; ++i) {
        x[i] = b[i];
        for (int j = A.GetRow()[i]; j < A.GetRow()[i + 1]; ++j) {
            x[i] -= A.GetVal()[j] * x[A.GetCol()[j]];
        }
        x[i] /= diag[i];
    	}
	}

    return x;
}


#endif