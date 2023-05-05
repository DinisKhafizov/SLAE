#ifndef SYM_GS_METH_WITH_CHEB_ACCEL
#define SYM_GS_METH_WITH_CHEB_ACCEL

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "Acceleration/Cheb_Accel.hpp"
#include "GS_Iterations.hpp"

std::vector<double> SGSMwCA(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_max);

#endif