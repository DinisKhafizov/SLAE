#ifndef PROB1
#define PROB1

#include <iostream>
#include <vector>
#include <string>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"

std::pair<double, double> lambda(const double a, const double b, const int L) {
    std::pair<double, double> res;
    res.first = 2 * (b - 2 * std::abs(a) * cos(M_PI/(L + 1)));
    res.second = 2 * (b + 2 * std::abs(a) * cos(M_PI/(L + 1)));
    return res;
}

std::vector<double> SIM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double a, const double B, const int L, std::string excersise) {
    double tau;
    std::pair<double, double> lam = lambda(a, B, L);
    if (excersise == "first") {
        tau = 1/lam.second;
    }
    else{
        tau = 2/(lam.first + lam.second);
    }
    std::vector<double> res;
    int count = 0;
	std::vector<double> x = x_0, vec(size(b)), r(size(b));

	while(first_norm(r) > tolerance) {
		vec = A * x;
        r = vec - b;
        count++;
        res.push_back(count);
        res.push_back(log(first_norm(r)));
        x = x - tau*r;
	}
	return res;
}

std::vector<double> SIMwCA(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, 
const double a, const double B, const int L, const int n) {

    double lambda_min, lambda_max;
    std::vector<double> res;
    std::pair<double, double> lam = lambda(a, B, L);
    lambda_min = lam.first;
    lambda_max = lam.second;
	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	int N = size(x_0); 
	double norm = 0; 
	std::vector<double> x = x_0;
    std::vector<double> vec = A * x;

	for (int i = 0; i < pow(2, n); ++i) {
		x = x - taus[nums[i]] * (vec - b);
        vec = A * x;
		norm = first_norm(vec - b);
	}

	if (norm > tolerance) {
		return SIMwCA(A, x, b, tolerance, a, B, L, n);
	}
	else {
		return x;
	}

}



#endif