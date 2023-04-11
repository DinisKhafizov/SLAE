#ifndef STEEPEST_GD
#define STEEPEST_GD

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SGD(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
	double tau;
	std::vector<double> x = x_0, r(size(b), 100);

	while(first_norm(r) > tolerance) {
		r = A*x - b;
		tau = (r * r) / (A * r * r);
		x = x - tau * r;
	}
	return x;
	
}

#endif