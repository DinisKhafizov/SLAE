#include "BiCG.hpp"

std::vector<double> BiCG(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, const double tolerance) {
    std::vector<double> r = A * x_0 - b;
    std::vector<double> r_t = r, p = r, p_t = r, x(A.GetN());
    double rho = r * r_t, rho_last;
    CSR A_t = A.transpose();
    std::vector<double> z = A * p, z_t = A_t * p_t;
    double q = rho/(p_t * z), theta;
    x = x_0 - q * p;
    r = r - q * z;
    if (first_norm(r) < tolerance) {
        return x;
    }
    r_t = r_t - q * z_t;
    rho_last = rho;

    while (first_norm(r) > tolerance) {
        rho = r_t * r;
        if (rho == 0) {
            return BiCG(A, b, x, tolerance);
        }
        theta = rho/rho_last;
        p = r + theta*p;
        p_t = r_t + theta*p_t;
        z = A*p;
        z_t = A_t * p_t;
        q = rho/(p_t * z);
        x = x - q*p;
        r = r - q*z;
        r_t = r_t - q*z_t;
    }
    return x;

}