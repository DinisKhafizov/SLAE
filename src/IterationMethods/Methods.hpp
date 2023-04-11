#ifndef ITER_METHODS
#define ITER_METHODS

#include <iostream>
#include <cmath>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Acceleration/Cheb_Accel.hpp"

std::vector<double> JM(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {

    int N = A.GetN(); //shape of matrix

    std::vector<double> x = x_0, diag(N), vec(N), x_last; // x and x_last for counting norm, vec for storage A*x, diag for storage diag els of A

    double norm = 100; // norm for comparison with tolerance that will be counted on each iteration

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
				A.GetVal()[k] = 0; 
			}
		}
	}
	

    while (norm > tolerance) {

		x_last = x;

        vec = A*x;

        for (int i = 0; i < N; ++i) {
            x[i] = (b[i] - vec[i]) / diag[i]; //new x
        }

		norm = 0;
		for (int i = 0; i < N; ++i) {
			norm += std::abs(x[i] - x_last[i]); //counting norm
		}

	
    }

    return x;


}

std::vector<double> GSM(CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {

	int N = A.GetN(), u; 

	double norm = 100, Ux = 0, Lx = 0; 
	
	std::vector<double> x = x_0, vec(N), x_last, diag(N);

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
			}
		}
	}

	while (norm > tolerance) {

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

		norm = 0;

        for (int i = 0; i < N; ++i) {
            norm += std::abs(x[i] - x_last[i]);
        }
	}

    return x;
}

std::vector<double> SIM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau) {

	int N = A.GetN();

	double norm = 0;
	
	std::vector<double> x = x_0, vec(N);

	while(norm > tolerance) {

		vec = A * x;

		for (int k = 0; k < N; ++k) {
			x[k] = x[k] - tau * (vec[k] - b[k]);
		}

		norm = 0;

        for (int i = 0; i < N; ++i) {
            norm += std::abs(vec[i] - b[i]);
		}

	}

	return x;
}

std::vector<double> SIMwCA(const CSR &A, const double lambda_min, const double lambda_max, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const int n) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);

	int N = A.GetN(); //Simple.

	double norm = 0; 
	
	std::vector<double> x = x_0, vec(N);

	for (int i = 0; i < pow(2, n); ++i) {

		vec = A * x;

		for (int k = 0; k < N; ++k) {
			x[k] = x[k] - taus[nums[i]] * (vec[k] - b[k]);
		}

		norm = 0;

        for (int i = 0; i < N; ++i) {
            norm += std::abs(vec[i] - b[i]);
		}

	}

	if (norm > tolerance) {
		return SIMwCA(A, lambda_min, lambda_max, x, b, tolerance, n);
	}
	else {
		return x;
	}

}

#endif