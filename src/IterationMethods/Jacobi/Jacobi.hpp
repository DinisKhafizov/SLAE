#ifndef JACOBI_METH
#define JACOBI_METH

#include <iostream>
#include <vector>
#include "CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> JM(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance) {

    const int N = size(x_0); //shape of matrix
    std::vector<double> x = x_0, diag(N), x_last(N, 100); //, vec(N); x and x_last for counting norm, vec for storage A*x, diag for storage diag els of A

	for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
				A.GetVal()[k] = 0; 
			}
		}
	}
    while (first_norm(x - x_last) > tolerance) {
		x_last = x;
        /*
        for (int i = 0; i < N; ++i) {
            x[i] = (b[i] - vec[i]) / diag[i]; //new x
        }
        */
        x = (b - A * x)/diag; //el-wise division
    }

    return x;
    
}

#endif