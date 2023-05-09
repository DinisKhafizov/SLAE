#include "Householder.hpp"

std::pair<Matrix, Matrix> QR_Decomposition(Matrix A) {
    const int N = A.GetN(), M = A.GetM();
    std::vector<double> v(M), nulls(M-1, 0);
    std::pair<Matrix, Matrix> QR;
    double koef;
    Matrix Q;

    v = A.getCol(0);
    if (v[0] > 0) { v[0] += second_norm(v); }
    else          { v[0] -= second_norm(v); }
    koef = -2/(v * v);
    Q = dot_self(v) ;
    Q * koef;
    Q.add_Identity();

    A.get_vals()[0] += koef * (v * A.getCol(0)) * v[0];
    A.change_Col(0, 1, nulls);
    for (size_t i = 1; i < N; ++i) {
        A.change_Col(i, 0, A.getCol(i) + koef * (v * A.getCol(i)) * v);
    }
    
    for (size_t i = 1; i < N; ++i) {
        v.pop_back();
        nulls.pop_back();
        v = A.getCol(i, i);
        if (v[0] > 0) { v[0] += second_norm(v); }
        else          { v[0] -= second_norm(v); }
        koef = -2/(v * v);

        A.get_vals()[i * (N + 1)] += koef * (v * A.getCol(i, i)) * v[0];
        A.change_Col(i, i+1, nulls);
        for (size_t j = i + 1; j < N; ++j) {
            A.change_Col(j, i, A.getCol(j, i) + koef * (v * A.getCol(j, i)) * v);
        }
        for (size_t j = 0; j < M; ++j) {
            Q.change_Col(j, i, Q.getCol(j, i) + koef * (v * Q.getCol(j, i)) * v);
        }
    
    }
    Q.transpose();
    QR.first = Q;
    QR.second = A;
    return QR;
}