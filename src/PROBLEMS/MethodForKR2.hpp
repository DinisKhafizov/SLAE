#ifndef PROB1
#define PROB1

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"
#include "IterationMethods/Gauss-Seidel/GS_Iterations.hpp"



void Write(const std::vector<double> &res, std::string name) {
    std::ofstream path;
    path.open(name);
    for (int i = 0; i < int(size(res)/2); ++i) {
        path << res[i*2] << ";" << res[i*2 + 1] << "\n";
    }
}

void Write(const std::pair<std::vector<double>, std::vector<double>> &res, std::string name) {
    std::ofstream path;
    path.open(name);
    for (int i = 0; i < size(res.first); ++i) {
        path << res.first[i] << ";" << res.second[i] << "\n";
    }
}

std::pair<double, double> lambda(const double a, const double b, const int L) {
    std::pair<double, double> res;
    res.first = 2 * (b - 2 * std::abs(a) * cos(M_PI/(L + 1)));
    res.second = 2 * (b + 2 * std::abs(a) * cos(M_PI/(L + 1)));
    return res;
}

/*
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
const double lambda_max, const int n, std::vector<double> &res, double counter = 0) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	const int N = size(x_0), it = pow(2, n);
    int i = 0;  
	std::vector<double> x = x_0;
    std::vector<double> r = A * x - b;
    double norm = first_norm(r);// = first_norm(r);
    //std::vector<double> res = {counter, log(norm)};

	//for (int i = 0; i < pow(2, n); ++i) {
    while ((i < it) && (norm > tolerance)) {
		x = x - taus[nums[i]] * r;      
        r = A * x - b;
		norm = first_norm(r);
        i++;
        counter += 1;
        res.push_back(counter);
        res.push_back(log(norm));
	}

	if (norm > tolerance) {
		return SIMwCA(A, x, b, tolerance, lambda_min, lambda_max, n, res, counter);
	}
	else {
		return res;
	}

}

std::vector<double> SIMwCA_next(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_min, 
const double lambda_max, const int n, double counter = 0) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	const int N = size(x_0), it = pow(2, n);
    int i = 0;  
	std::vector<double> x = x_0;
    std::vector<double> r = A * x - b, res;
    double norm = first_norm(r);// = first_norm(r);
    //std::vector<double> res = {counter, log(norm)};

	//for (int i = 0; i < pow(2, n); ++i) {
    while ((i < it) && (norm > tolerance)) {
		x = x - taus[nums[i]] * r;      
        r = A * x - b;
		norm = first_norm(r);
        i++;
        counter += 1;
	}

	if (norm > tolerance) {
		return SIMwCA_next(A, x, b, tolerance, lambda_min, lambda_max, n, counter);
	}
	else {
        res.push_back(counter);
        res.push_back(lambda_max);
		return res;
	}

}



std::vector<double> SGSM(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {
	const int N = size(x_0); 
	std::vector<double> x = x_0, x_last(N, 100), diag(N);
    std::vector<double> r = A*x - b;
    double counter = 0;
    std::vector<double> res = {counter, log(first_norm(r))};

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
				A.GetVal()[k] = 0;
			}
		}
	}
	while(first_norm(r) > tolerance) {
		x_last = x;
		x = TopDownIteration(A, diag, b, x);
		x = DownUpIteration(A, diag, b, x);
        counter += 1;
        r = A*x - b + elWise_Mult(diag, x);
        res.push_back(counter);
        res.push_back(log(first_norm(r)));
	}
    return res;
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
    //r = A * x - b + elWiseMult(diag, x);
    //norm = first_norm(r);
    //counter += 1;
    //res.push_back(counter);
    //res.push_back(norm);

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
        r = A * x - b + elWise_Mult(diag, x);
        norm = first_norm(r);
        counter += 1;
        res.push_back(counter);
        res.push_back(log(norm));

    }
    return res;
}
*/
std::pair<std::vector<double>, std::vector<double>> SIM_in_proj(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau) {

	std::vector<double> x = x_0;
    std::vector<double> r = A * x - b;
    double norm = first_norm(r);
    std::pair<std::vector<double>, std::vector<double>> res;
    res.first.push_back(x_0[0]); 
    res.second.push_back(x_0.back());

	while(norm > tolerance) {
        x = x - tau*r;
        r = A * x - b;
        norm = first_norm(r);
        res.first.push_back(x[0]); 
        res.second.push_back(x.back());
	}

	return res;
}

std::pair<std::vector<double>, std::vector<double>> SGD_in_proj(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
	double tau;
	std::vector<double> x = x_0, r = A * x_0 - b;
    std::pair<std::vector<double>, std::vector<double>> res;
    res.first.push_back(x_0[0]); 
    res.second.push_back(x_0.back());

	while(first_norm(r) > tolerance) {
		tau = (r * r) / (A * r * r);
		x = x - tau * r;
		r = A*x - b;
        res.first.push_back(x[0]);
        res.second.push_back(x.back());
	}
	return res;	
}
std::vector<double> SIMwCA_in_proj(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_min, 
const double lambda_max, const int n, std::vector<double> res = {0}) {

	std::vector<int> nums = Numbers(n);
	std::vector<double> taus = TauFromCheb(n, lambda_min, lambda_max);
	const int N = size(x_0), it = pow(2, n);
	int i = 0; 
	double norm = first_norm(A * x_0 - b); 
	std::vector<double> x = x_0;
    std::vector<double> r = A*x - b;
    if (size(res) == 1) {
        res[0] = x_0[0];
        res[1] = x_0.back();
    }

	while ((i < it) && (norm > tolerance)) {
		x = x - taus[nums[i]] * r;
        r = A * x - b;
		norm = first_norm(r);
		++i;
        res.push_back(x[0]);
        res.push_back(x.back());
	}

	if (norm > tolerance) {
		return SIMwCA_in_proj(A, x, b, tolerance, lambda_min, lambda_max, n, res);
	}
	else {
		return res;
	}

}

std::pair<std::vector<double>, std::vector<double>> CGD_in_proj(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
    std::vector<double> x = x_0, r = A * x_0 - b;
    std::vector<double> d = r;
    double alpha, koef1 = r * r, koef2;
    std::pair<std::vector<double>, std::vector<double>> res;
    res.first.push_back(x_0[0]); 
    res.second.push_back(x_0.back());
    while(first_norm(r) > tolerance) {
        alpha = (r * r)/(d * (A * d));
        x = x - alpha * d;
        r = A * x - b;
        koef2 = r * r;
        d = r + d * koef2/koef1;
        koef1 = koef2;
        res.first.push_back(x[0]);
        res.second.push_back(x.back()); 
        }
    return res;
}






#endif