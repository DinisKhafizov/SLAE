#include "SGD.hpp"

std::vector<double> SGD(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
	double tau;
	std::vector<double> x = x_0, r = A * x_0 - b;

	while(first_norm(r) > tolerance) {
		tau = (r * r) / (A * r * r);
		x = x - tau * r;
		r = A*x - b;
	}
	return x;
	
}