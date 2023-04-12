#ifndef PROB1
#define PROB1

#include <iostream>
#include <vector>
#include <fstream>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"
#include "IterationMethods/Gauss-Seidel/GS_Iterations.hpp"

void Write(const std::vector<double> &res) {
    std::ofstream path;
    path.open("Problem1.txt");
    for (int i = 0; i < int(size(res)/2); ++i) {
        path << res[i*2] << ";" << res[i*2 + 1] << "\n";
    }
}

std::pair<double, double> lambda(const double a, const double b, const int L) {
    std::pair<double, double> res;
    res.first = 2 * (b - 2 * std::abs(a) * cos(M_PI/(L + 1)));
    res.second = 2 * (b + 2 * std::abs(a) * cos(M_PI/(L + 1)));
    return res;
}

std::vector<double> SIM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau) {

	std::vector<double> x = x_0;
    std::vector<double> r = A * x - b;
    double norm = first_norm(r), counter = 0;
    std::vector<double> res = {counter, log(norm)};

	while(norm > tolerance) {
        x = x - tau*r;
        r = A * x - b;
        norm = first_norm(r);
        counter += 1;
        res.push_back(counter);
        res.push_back(log(norm));
	}

	return res;
}

std::vector<double> SIMwCA(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_min, 
const double lambda_max, const int n, double counter = 0) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	int N = size(x_0);  
	std::vector<double> x = x_0;
    std::vector<double> r = A * x - b;
    double norm = first_norm(r);
    std::vector<double> res = {counter, log(norm)};

	for (int i = 0; i < pow(2, n); ++i) {
		x = x - taus[nums[i]] * r;
        r = A * x - b;
		norm = first_norm(r);
        counter += 1;
        res.push_back(counter);
        res.push_back(log(norm));
	}

	if (norm > tolerance) {
		return SIMwCA(A, x, b, tolerance, lambda_min, lambda_max, n, counter);
	}
	else {
		return res;
	}

}

std::vector<double> SGSMwCA(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_max) {
    const int N = size(x_0); 
	double mu_0 = 0, mu_1 = 1/lambda_max, mu_2; 
	std::vector<double> x = x_0, x_last, diag(N), x_last2;
    std::vector<double> r = A * x - b;
    double norm = first_norm(r), counter = 0;
    std::vector<double> res = {counter, log(norm)};

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
    norm = first_norm(r);
    counter += 1;
    res.push_back(counter);
    res.push_back(norm);

    while (norm > tolerance){  
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
        norm = first_norm(r);
        counter += 1;
        res.push_back(counter);
        res.push_back(log(norm));

    }
    return res;
}




#endif