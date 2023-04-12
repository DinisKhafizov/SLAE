#ifndef THREE_DIAG_SOLVER_TRIDIAG_SOLVER_H
#define THREE_DIAG_SOLVER_TRIDIAG_SOLVER_H


#include <iostream>
#include <vector>
#include <string>
#include <fstream>


class TridiagonalMatrix 
{
private:
    int n;
    std::vector<double> A, B, C;
public:
    TridiagonalMatrix(const std::vector<double> &A1, const std::vector<double> &B1,const std::vector<double> &C1): A{A1}, B{B1}, C{C1}, n{size(B1)}{}

    int GetN() {
        return n;
    }
    double GetA(int i) {
        return A[i];
    }
    double GetB(int i) {
        return B[i];
    }
    double GetC(int i) {
        return C[i];
    }
};

#endif //THREE_DIAG_SOLVER_TRIDIAG_SOLVER_H