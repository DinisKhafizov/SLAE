#ifndef BICG
#define BICG

#include "Vect/VectorOperations.hpp"

std::vector<double> BiCG(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, const double tolerance);

#endif