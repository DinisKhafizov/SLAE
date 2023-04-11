#ifndef STEEPEST_GD
#define STEEPEST_GD

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SGD(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, double tolerance) {
	double tau;
	std::vector<double> x = x_0, vec(size(b)), r(size(b));

	while(first_norm(r) > tolerance) {
		vec = A * x;
		r = vec - b;
		tau = (r * r) / (A * r * r);
		x = x - tau * r;
	}
	return x;
	
}

#endif