#ifndef SOR_MY
#define SOR_MY

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SOR(CSR A, const std::vector<double> x_0, std::vector<double> &b, const double tolerance, const double omega) {

	int N = A.GetN(), u; 

	double Ux = 0, Lx = 0, x_i, koef = 1 - omega; 
	
	std::vector<double> x = x_0, vec(N), x_last(N, 100), diag(N);

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
			x_i = x[i];
			x[i] = omega * b[i];
			for (int j = A.GetRow()[i]; j < A.GetRow()[i + 1]; ++j) {
				x[i] -= omega * A.GetVal()[j] * x[A.GetCol()[j]];
			}
			x[i] /= diag[i];
			x[i] += koef * x_i;
    	}
			//x[j] = (1 - omega) * x[j] + omega * (b[j] - Ux - Lx) / diag[j];
	}

	return x;

}



#endif