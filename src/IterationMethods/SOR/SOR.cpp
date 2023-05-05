#include "SOR.hpp"

std::vector<double> SOR(const CSR &A, const std::vector<double> &x_0, std::vector<double> b, const double tolerance, const double omega) {
	const int N = size(x_0);
	int counter, row_begin, row_end;
	double x_i, koef = 1 - omega;
	std::vector<double> x = x_0, x_last(N, 100), diag(N);

	diag = A.getDiag();
	b = omega * b;

	while (first_norm(x - x_last) > tolerance) {
		x_last = x;
		for (size_t i = 0; i < N; ++i) {
			x_i = x[i];
        	x[i] = b[i];
        	counter = 1;
        	row_begin = A.GetRow()[i];
        	row_end = A.GetRow()[i + 1];
        	for (size_t j = row_begin; j < row_end && A.GetCol()[j] < i; ++j) {
            	x[i] -= omega * A.GetVal()[j] * x[A.GetCol()[j]];
            	++counter;
        	}
        	for (size_t j = row_begin + counter; j < row_end; ++j) {
            	x[i] -= omega * A.GetVal()[j] * x[A.GetCol()[j]];
        	} 
        	x[i] /= diag[i];
			x[i] += koef * x_i;
    	}
	}

	return x;

}