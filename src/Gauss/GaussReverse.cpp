#include "GaussReverse.hpp"

std::vector<double> GaussReverse(const Matrix &A, std::vector<double> b) {
    const int N = A.GetN(); 
    const int N_pOne = N + 1;
    int Ni = N, Nvals = N * N- 1;
    double sub;
    for (int i = 0, end=N-1; i < end; ++i) {
        Ni -= 1; 
        b[Ni] /= A.get_vals()[Nvals];
 
        for (int j = 1, end = N-i; j < end; ++j) {
            sub = A.get_vals()[Nvals - N * j];
            b[Ni - j] -= sub * b[Ni];
        }
        Nvals -= N_pOne;
    }
    b[0] /= A.get_vals()[0];
    return b;
}
std::vector<double> GaussReverse(const Matrix &A, std::vector<double> b, const int minor) {
    const int N = A.GetN();
    const int N_pOne = N + 1;
    int Ni = minor, Nvals = minor + N * (minor - 1) - 1;
    double sub;
    for (int i = 0, end=minor-1; i < end; ++i) {
        Ni -= 1; 
        b[Ni] /= A.get_vals()[Nvals];
        A.get_vals()[Nvals] = 1;
        for (int j = 1, end = minor-i; j < end; ++j) {
            sub = A.get_vals()[Nvals - N * j];
            b[Ni - j] -= sub * b[Ni];
        }
        Nvals -= N_pOne;
    }
    b[0] /= A.get_vals()[0];
    return b;
}