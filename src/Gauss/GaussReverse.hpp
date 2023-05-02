#ifndef GAUSS_REV
#define GAUSS_REV

#include "Matrixes/SMatrix/StandardMatrix.hpp"

std::vector<double> GaussReverse(Matrix A, std::vector<double> b) {
    const int N = A.GetN(); // N_mOne = A.GetN(), N_pOne = A.GetN() + 1, Nsqr = A.GetN() * A.GetN(), Nsqr_mOne = A.GetN() * A.GetN() - 1;
    const int N_pOne = N + 1, Nsqr = N*N;
    int Ni = N, Nvals = Nsqr - 1;
    double sub;
    for (int i = 0, end=N-1; i < end; ++i) {
        Ni -= 1; 
        b[Ni] /= A.get_vals()[Nvals];
        A.get_vals()[Nvals] = 1;
        for (int j = 1, end = N-i; j < end; ++j) {
            sub = A.get_vals()[Nvals - N * j];
            A.get_vals()[Nvals - N*j] = 0;
            b[Ni - j] -= sub * b[Ni];
        }
        Nvals -= N_pOne;
    }
    b[0] /= A.get_vals()[0];
    return b;
}

#endif