#ifndef ITER_METHODS_PROBLEM3
#define ITER_METHODS_PROBLEM3

#include <iostream>
#include <cmath>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include <fstream>

/*Difference between this file and other methods is in: 
Methods' not returning x, it returns only Number of iterations.*/

int SimpleIterationMethod(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double tau){

	int N = A.GetN(), num_iter = 0; //Simple.

	double norm = 100, Ux = 0, Lx = 0; 
	
	std::vector<double> x = x_0, vec(N);

	while(norm > tolerance) {


		vec = A * x;

		for (int k = 0; k < N; ++k) {
			x[k] = x[k] - tau * (vec[k] - b[k]);
		}

		norm = 0;

        for (int i = 0; i < N; ++i) {
            norm += std::abs(vec[i] - b[i]);
		}

		num_iter++;
}

	return num_iter;
}


#endif