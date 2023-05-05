#ifndef TEST_THREE_DIAG_SOLVER_H
#define TEST_THREE_DIAG_SOLVER_H

#include <iostream>
#include <vector>
#include "Matrixes/Tridiag/Tridiag.hpp"




std::vector<double> Progonka(TridiagonalMatrix Matrix, std::vector<double> D);

#endif