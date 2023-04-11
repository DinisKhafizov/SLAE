#ifndef CHEBISHEV_ACCELERATION
#define CHEBISHEV_ACCELERATION

#include <iostream>
#include <vector>
#include <cmath>

std::vector<int> Numbers(int n) {
    if (n == 1) {
        std::vector<int> a = {0, 1};
        return a;
    }
    int two_pow_n = pow(2, n), two_pow_i = 1;
    std::vector<int> numb(two_pow_n), prev(two_pow_n);
    prev[0] = 0;
    prev[1] = 1;
    for (int i = 1; i < n; ++i) {
        two_pow_i *= 2;
        for (int j = 0; j < two_pow_i; ++j) {
            numb[j*2] = prev[j];
            numb[j*2 + 1] = two_pow_i * 2 - 1 - prev[j];
        }
        for (int j = 0; j < two_pow_i*2; ++j) {
            prev[j] = numb[j];
        }
    }
    return numb;
}

std::vector<double> TauFromCheb(const int n, const double lambda_min, const double lambda_max) {
    std::vector<double> taus(n);
    const double cos_const = cos(M_PI/n), sin_const = sin(M_PI/n), lam_plus = (lambda_max + lambda_min)/2, lam_minus = (lambda_max - lambda_min)/2;
    double Sin = sin(M_PI/(2 * n));
    taus[0] = cos(M_PI/(2 * n));
    for(int i = 1; i < n; ++i) {
        taus[i] = taus[i - 1] * cos_const - Sin * sin_const;
        Sin = Sin * cos_const + taus[i - 1] * sin_const;
    }
    for (int i = 0; i < n; ++i) {
        taus[i] = 1 / (lam_plus + lam_minus * taus[i]);
    }
    return taus;
}

#endif