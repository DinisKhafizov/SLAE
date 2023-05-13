#ifndef CONJUGATE_GRADIENT_METH
#define CONJUGATE_GRADIENT_METH

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> CGD(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance);
std::vector<double> CGD_precond(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance);

#endif