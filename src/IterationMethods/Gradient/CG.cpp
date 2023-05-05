#include "CG.hpp"

std::vector<double> CGD(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
    std::vector<double> x = x_0, r = A * x_0 - b;
    std::vector<double> d = r;
    double alpha, koef1 = r * r, koef2;
    while(first_norm(r) > tolerance) {
        alpha = (d * r)/(d * (A * d));
        x = x - alpha * d;
        r = A * x - b;
        koef2 = r * r;
        d = r + d * koef2/koef1;
        koef1 = koef2; 
        }
    return x;
}