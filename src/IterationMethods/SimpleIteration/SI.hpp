#ifndef SIMPLE_ITERATION_METH
#define SIMPLE_ITERATION_METH

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SIM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau);

#endif