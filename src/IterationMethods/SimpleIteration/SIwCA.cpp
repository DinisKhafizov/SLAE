#include "SIwCA.hpp"

std::vector<double> SIMwCA(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_min, 
const double lambda_max, const int n) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	const int N = size(x_0), it = pow(2, n);
	int i = 0; 
	double norm = 0; 
	std::vector<double> x = x_0;
    std::vector<double> r = A*x - b;

	while (i < it || norm > tolerance) {
		x = x - taus[nums[i]] * r;
        r = A * x - b;
		norm = first_norm(r);
		++i;
	}

	if (norm > tolerance) {
		return SIMwCA(A, x, b, tolerance, lambda_min, lambda_max, n);
	}
	else {
		return x;
	}

}