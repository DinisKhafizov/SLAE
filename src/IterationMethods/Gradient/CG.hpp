#ifndef CONJUGATE_GRADIENT_METH
#define CONJUGATE_GRADIENT_METH

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> CGD(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
    std::vector<double> x = x_0, r = A * x_0 - b;
    std::vector<double> d = r, r_last = r;
    double alpha, koef1 = r * r, koef2;
    while(first_norm(r) > tolerance) {
        alpha = (r * r)/(d * (A * d));
        x = x - alpha * d;
        r = A * x - b;
        koef2 = r * r;
        d = r + d * koef2/koef1;
        koef1 = koef2; 
        }
    return x;
}

#endif