#ifndef CHEBISHEV_ACCELERATION
#define CHEBISHEV_ACCELERATION

#include <iostream>
#include <vector>
#include <cmath>

std::vector<int> Numbers(const int n);

std::vector<double> TauFromCheb(const int n, const double lambda_min, const double lambda_max);

#endif