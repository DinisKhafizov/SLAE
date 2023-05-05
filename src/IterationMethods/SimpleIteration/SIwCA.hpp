#ifndef SIMPLE_ITERATION_METH_WITH_CHEB_ACCEL
#define SIMPLE_ITERATION_METH_WITH_CHEB_ACCEL

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"

std::vector<double> SIMwCA(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_min, 
const double lambda_max, const int n);
#endif