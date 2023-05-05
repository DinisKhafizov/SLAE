#include "SGSwCA.hpp"

std::vector<double> SGSMwCA(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance, const double lambda_max) {
    const int N = size(x_0); 
	double mu_0 = 0, mu_1 = 1/lambda_max, mu_2; 
	std::vector<double> x = x_0, x_last, diag(N), x_last2;
    std::vector<double> r;

    diag = A.getDiag();

    x_last = x;
    x = TopDownIteration(A, diag, b, x);
    x = DownUpIteration(A, diag, b, x);
    r = A * x - b;

    while (first_norm(x - x_last) > tolerance){  
        x_last2 = x_last;
        x_last = x;
        x = TopDownIteration(A, diag, b, x);
        x = DownUpIteration(A, diag, b, x);
        mu_2 = 2 * mu_1 / lambda_max - mu_0;

        x = (2*mu_1/lambda_max/mu_2) * x;
        x = x - (mu_0/mu_2) * x_last2;
        mu_0 = mu_1;
        mu_1 = mu_2; 
        r = A * x - b;
    }
    return x;
}