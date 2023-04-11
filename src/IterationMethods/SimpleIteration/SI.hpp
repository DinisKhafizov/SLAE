#ifndef SIMPLE_ITERATION_METH
#define SIMPLE_ITERATION_METH

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SIM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau) {

	std::vector<double> x = x_0, vec(size(b));

	while(first_norm(vec - b) > tolerance) {
		vec = A * x;
        /*
		for (int k = 0; k < N; ++k) {
			x[k] = x[k] - tau * (vec[k] - b[k]);
		}
        */
        x = x - tau*(vec - b);
	}

	return x;
}

#endif