//#ifndef QR_MATRIX
//#define QR_MATRIX

#include <iostream>
#include <vector>
#include "Matrixes/SMatrix/StandardMatrix.hpp"
#include "Vect/VectorOperations.hpp"

std::pair<Matrix, Matrix> QR_Decomposition(const Matrix &A) {
    const int N = A.GetN();
    double koef;
    Matrix Q, R;
    std::vector<double> v;

    v = A.getCol(0);
    if (v[0] > 0) {
        v[0] += second_norm(v);
    }
    else {
        v[0] -= second_norm(v);
    }
    koef = -2/(v * v);
    Q = dot_self(v);
    Q * koef;


    for (size_t i = 1; i < N; ++i) {
        v = A.getCol(i);
        if (v[0] > 0) {
            v[0] += second_norm(v);
        }
        else {
            v[0] -= second_norm(v);
        }
        koef = 2/(v * v);



    }
}

//#endif 
