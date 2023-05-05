#ifndef SOR_MY
#define SOR_MY

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

std::vector<double> SOR(const CSR &A, const std::vector<double> &x_0, std::vector<double> b, const double tolerance, const double omega);


#endif