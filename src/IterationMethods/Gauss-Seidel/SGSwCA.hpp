#ifndef SYM_GS_METH_WITH_CHEB_ACCEL
#define SYM_GS_METH_WITH_CHEB_ACCEL

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"
#include "GS_Iterations.hpp"

std::vector<double> SGSMwCA(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_max) {
    const int N = size(x_0); 
	double mu_0 = 0, mu_1 = 1/lambda_max, mu_2; 
	std::vector<double> x = x_0, x_last, diag(N), x_last2;
    std::vector<double> r;

    for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { 
				diag[i] = A.GetVal()[k];
                A.GetVal()[k] = 0;
			}
		}
	}

    x_last = x;
    x = TopDownIteration(A, diag, b, x);
    x = DownUpIteration(A, diag, b, x);
    r = A * x - b;

    while (first_norm(x - x_last) > tolerance){  
        x_last2 = x_last;
        x_last = x;
        x = TopDownIteration(A, diag, b, x);
        x = DownUpIteration(A, diag, b, x);
        mu_2 = 2 * mu_1 / lambda_max - mu_0;

        x = (2*mu_1/lambda_max/mu_2) * x;
        x = x - (mu_0/mu_2) * x_last2;
        mu_0 = mu_1;
        mu_1 = mu_2; 
        r = A * x - b;
    }
    return x;
}

#endif