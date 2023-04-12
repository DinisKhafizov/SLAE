#ifndef SYM_GS_METH
#define SYM_GS_METH

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "GS_Iterations.hpp"

std::vector<double> SGSM(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {
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
	while(first_norm(x - x_last) > tolerance) {
		x_last = x;
		x = TopDownIteration(A, diag, b, x);
		x = DownUpIteration(A, diag, b, x);
	}
    return x;
}
#endif