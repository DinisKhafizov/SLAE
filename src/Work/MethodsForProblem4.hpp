#ifndef ITER_METHODS_PROBLEM4
#define ITER_METHODS_PROBLEM4

#include <iostream>
#include <cmath>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include <fstream>


/*Difference between this file and other methods is in: 
Methods return not x, but pair of vectors: 1-st is vector of Number of iterations and 2-nd is vector of norm between (b - A*x)
*/
	

std::pair<std::vector<double>, std::vector<int>> JacobiMethod(CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {

    int N = A.GetN(), y = 0; //shape of matrix

    std::vector<double> x = x_0, diag(N), vec(N), r; 
	
	std::vector<int> num_iter;

    double norm = 100; // norm for comparison with tolerance that will be counted on each iteration

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
				A.GetVal()[k] = 0; 
			}
		}
	}
	

    while (norm > tolerance) {

        vec = A*x;

        norm = 0;

		for (int i = 0; i < N; ++i) {
            norm += std::abs(vec[i] + diag[i] * x[i] - b[i]);
        }

        for (int i = 0; i < N; ++i) {
            vec[i] = (b[i] - vec[i]) / diag[i];
        }

		x = vec;


		y++;
		num_iter.resize(y);
		num_iter[y - 1] = y;
		r.resize(y);
		r[y - 1] = norm;
    }

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) {
				A.GetVal()[k] = diag[i];
			}
		}
	}//returning A its values

	std::pair<std::vector<double>, std::vector<int>> l(r, num_iter);

    return l;


}

std::pair<std::vector<double>, std::vector<int>> GaussSeidelMethod(CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {

	int N = A.GetN(), u, y = 0;

	double norm = 100, Ux = 0, Lx = 0;
	
	std::vector<double> x = x_0, vec(N), r;
	std::vector<int> num_iter;

	while (norm > tolerance) {

		for (int k = 1; k < A.GetRow()[1]; ++k) {
			Ux += A.GetVal()[k] *  x[A.GetCol()[k]];
		}

		x[0] = (b[0] - Ux) / A(0, 0);

		for (int j = 1; j < N; ++j) {

			Ux = 0; //sum U * x_i
			Lx = 0; // sum L * x_(i+1)
			u = 1;

			for (int k = A.GetRow()[j]; k < A.GetRow()[j + 1] && A.GetCol()[k] < j; ++k) {
				Lx += A.GetVal()[k] * x[A.GetCol()[k]];
				u++; //counter of iterations (A - square matrix, so *num of counting Ux iterations*  = N - *number of counting Lx iterations*)
			}

			for (int k = A.GetRow()[j] + u; k < A.GetRow()[j + 1]; ++k) { //so here I use counter u
				Ux += A.GetVal()[k] * x[A.GetCol()[k]];
			}

			x[j] = (b[j] - Ux - Lx) / A(j, j);

		}

		vec = A * x;

		norm = 0;

        for (int i = 0; i < N; ++i) {
            norm += std::abs(vec[i] - b[i]);
        }
		y++;
		num_iter.resize(y);
		num_iter[y - 1] = y;
		r.resize(y);
		r[y - 1] = norm;
	}
	std::pair<std::vector<double>, std::vector<int>> l(r, num_iter);

    return l;
}

std::pair<std::vector<double>, std::vector<int>> SimpleIterationMethod(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double tau){

	int N = A.GetN(), y = 0; //Simple.

	double norm = 100, Ux = 0, Lx = 0; 
	
	std::vector<double> x = x_0, vec(N), r;

	std::vector<int> num_iter;

	while(norm > tolerance) {


		vec = A * x;

		for (int k = 0; k < N; ++k) {
			x[k] = x[k] - tau * (vec[k] - b[k]);
		}

		norm = 0;

        for (int i = 0; i < N; ++i) {
            norm += std::abs(vec[i] - b[i]);
		}

		y++;
		num_iter.resize(y);
		num_iter[y - 1] = y;
		r.resize(y);
		r[y - 1] = norm;
	
}
	std::pair<std::vector<double>, std::vector<int>> l(r, num_iter);
	return l;
}

#endif //ITER_METHODS_PROBLEM4