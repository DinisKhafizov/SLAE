#ifndef GAUSS_REV
#define GAUSS_REV

#include "Matrixes/SMatrix/StandardMatrix.hpp"

std::vector<double> GaussReverse(const Matrix &A, std::vector<double> b);
std::vector<double> GaussReverse(const Matrix &A, std::vector<double> b, const int minor);

#endif