#ifndef SYM_GS_METH_WITH_CHEB_ACCEL
#define SYM_GS_METH_WITH_CHEB_ACCEL

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"
#include "GS_Iterations.hpp"

std::vector<double> SGSMwCA(CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_max) {
    int N = A.GetN(), u; 

	double Ux = 0, Lx = 0, mu_0 = 0, mu_1 = 1/lambda_max, mu_2; 
	
	std::vector<double> x = x_0, vec(N), x_last, diag(N), x_last2;

    for (int i = 0; i < N; ++i) {
		for(int k = A.GetRow()[i]; k < A.GetRow()[i + 1]; ++k) {
			if (A.GetCol()[k] == i) { //nulling diag elements of A (without changing Rows and Cols vectors of CSR) and writing them to vector diag
				diag[i] = A.GetVal()[k];
			}
		}
	}
    x_last = x;
    /*
    for (int k = 1; k < A.GetRow()[1]; ++k) {
        Ux += A.GetVal()[k] *  x[A.GetCol()[k]];
    }

    x[0] = (b[0] - Ux) / diag[0];

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

        x[j] = (b[j] - Ux - Lx) / diag[j];

    }

    Ux = 0;
    Lx = 0;

    for (int k = A.GetRow()[N - 1]; k < A.GetRow()[N] - 1; ++k) {
        Lx += A.GetVal()[k] * x[A.GetCol()[k]];
    }

    x.back() = (b[N - 1] - Lx) / diag[N - 1];

    for (int j = N - 1; j > 0; --j) {

        Ux = 0; //sum U * x_i
        Lx = 0; // sum L * x_(i+1)
        u = 1;

        for (int k = A.GetRow()[j - 1]; k < A.GetRow()[j] && A.GetCol()[k] < (j - 1); ++k) {
            Lx += A.GetVal()[k] * x[A.GetCol()[k]];
            u++; //counter of iterations (A - square matrix, so *num of counting Ux iterations*  = N - *number of counting Lx iterations*)
        }

        for (int k = A.GetRow()[j - 1] + u; k < A.GetRow()[j]; ++k) { //so here I use counter u
            Ux += A.GetVal()[k] * x[A.GetCol()[k]];
        }

        x[j - 1] = (b[j - 1] - Ux - Lx) / diag[j - 1];

    }
    */
    x = TopDownIteration(A, diag, b, x);
    x = DownUpIteration(A, diag, b, x);

    while (first_norm(x - x_last) > tolerance){
        
        x_last2 = x_last;
        x_last = x;
        /*
        Ux = 0; 
        Lx = 0;

        for (int k = 1; k < A.GetRow()[1]; ++k) {
            Ux += A.GetVal()[k] *  x[A.GetCol()[k]];
        }

        x[0] = (b[0] - Ux) / diag[0];

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

            x[j] = (b[j] - Ux - Lx) / diag[j];

        }

        Ux = 0;
        Lx = 0;

        for (int k = A.GetRow()[N - 1]; k < A.GetRow()[N] - 1; ++k) {
            Lx += A.GetVal()[k] * x[A.GetCol()[k]];
        }

        x.back() = (b[N - 1] - Lx) / diag[N - 1];

        for (int j = N - 1; j > 0; --j) {

            Ux = 0; //sum U * x_i
            Lx = 0; // sum L * x_(i+1)
            u = 1;

            for (int k = A.GetRow()[j - 1]; k < A.GetRow()[j] && A.GetCol()[k] < (j-1); ++k) {
                Lx += A.GetVal()[k] * x[A.GetCol()[k]];
                u++; //counter of iterations (A - square matrix, so *num of counting Ux iterations*  = N - *number of counting Lx iterations*)
            }
            for (int k = 1; k < A.GetRow()[1]; ++k) {
                Ux += A.GetVal()[k] *  x[A.GetCol()[k]];
            }

            x[j - 1] = (b[j - 1] - Ux - Lx) / diag[j - 1];
        }
        */

        x = TopDownIteration(A, diag, b, x);
        x = DownUpIteration(A, diag, b, x);
        mu_2 = 2 * mu_1 / lambda_max - mu_0;

        x = (2*mu_1/lambda_max/mu_2) * x;
        x = x - (mu_0/mu_2) * x_last2;
        mu_0 = mu_1;
        mu_1 = mu_2; 
    }

    return x;
}

#endif