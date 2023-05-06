#ifndef VECTOR_OP
#define VECTOR_OP

#include <iostream>
#include <vector>
#include <cmath>
#include "Matrixes/CSR/MatrixOnCSR.hpp"
#include "Matrixes/SMatrix/StandardMatrix.hpp"


/*
Three sections:
1) Arithmetic operations
2) Norms counting
3) Specific operations with vector

*/

// Arithmetic: '+', '-', '*', '/', '='

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> operator+(const std::vector<double> &a, double b);
std::vector<double> operator+=(const std::vector<double> &a, const std::vector<double> &b);

std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> operator-(const std::vector<double> &a, double b);
std::vector<double> operator-=(const std::vector<double> &a, const std::vector<double> &b);

double operator*(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> operator*(double x, const std::vector<double> &a);
std::vector<double> operator*(const std::vector<double> &a, double x);

std::vector<double> operator/(const std::vector<double> &a, double x);
std::vector<double> operator/(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> operator/=(const std::vector<double> &x, const double a);

bool operator==(const std::vector<double> &a, std::vector<double> &b);

// norms: first, second

double first_norm(const std::vector<double> &x);

double second_norm(const std::vector<double> &x);
double energetic_norm(const std::vector<double> &x, const CSR &A);

//specific

std::vector<double> elWise_Mult(const std::vector<double> &a, const std::vector<double> &b);
Matrix dot(const std::vector<double> &a, const std::vector<double> &b);
Matrix dot_self(const std::vector<double> &a);

#endif