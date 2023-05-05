#ifndef SYM_GS_METH
#define SYM_GS_METH

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "GS_Iterations.hpp"

std::vector<double> SGSM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance);
#endif