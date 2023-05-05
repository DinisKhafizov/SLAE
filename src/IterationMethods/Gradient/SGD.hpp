#ifndef STEEPEST_GD
#define STEEPEST_GD

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SGD(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance);

#endif