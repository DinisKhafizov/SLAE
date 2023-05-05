#ifndef JACOBI_METH
#define JACOBI_METH

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> JM(CSR A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance);

#endif