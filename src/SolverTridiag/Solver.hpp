#ifndef TEST_THREE_DIAG_SOLVER_H
#define TEST_THREE_DIAG_SOLVER_H

#include <iostream>
#include <vector>
#include "Matrixes/Tridiag/Tridiag.hpp"




std::vector<double> Progonka(TridiagonalMatrix Matrix, std::vector<double> D) 
{
    std::vector<double> x, p, q;
    double denom;
    const int N = Matrix.GetN();
    x.reserve(N);
    p.reserve(N - 1);
    q.reserve(N - 1);
    p[0] = - Matrix.GetC(0) / Matrix.GetB(0);
    q[0] = D[0] / Matrix.GetB(0);
    for (int i = 0; i < N - 2; ++i) {
        denom = Matrix.GetA(i) * p[i] + Matrix.GetB(i + 1);
        p[i + 1] = - Matrix.GetC(i + 1) / denom;
        q[i + 1] = (D[i + 1] - Matrix.GetA(i) * q[i]) / denom;
    }
    x[N - 1] = (D[N - 1] - Matrix.GetA(N - 2) * q[N - 2]) / (Matrix.GetA(N - 2) * p[N - 2] + Matrix.GetB(N - 1));
    for (int i = 1; i < N; ++i) {
        x[N - i - 1] = x[N - i] * p[N - i - 1] + q[N - i - 1];
    } 
    return(x);
}

#endif