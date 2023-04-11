#ifndef SIMPLE_ITERATION_METH_WITH_CHEB_ACCEL
#define SIMPLE_ITERATION_METH_WITH_CHEB_ACCEL

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"

std::vector<double> SIMwCA(const CSR &A, const double lambda_min, const double lambda_max, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const int n) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	int N = A.GetN(); 
	double norm = 0; 
	std::vector<double> x = x_0;
    std::vector<double> vec = A * x;

	for (int i = 0; i < pow(2, n); ++i) {

		for (int k = 0; k < N; ++k) {
			x[k] = x[k] - taus[nums[i]] * (vec[k] - b[k]);
		}

        vec = A * x;
		norm = first_norm(vec - b);
	}

	if (norm > tolerance) {
		return SIMwCA(A, lambda_min, lambda_max, x, b, tolerance, n);
	}
	else {
		return x;
	}

}

#endif