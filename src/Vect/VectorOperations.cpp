#include "VectorOperations.hpp"

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i){
        res[i] = a[i] + b[i];
    }
    return res;
}
std::vector<double> operator+(const std::vector<double> &a, double b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] + b;
    }
    return res;
}

std::vector<double> operator+=(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i){
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i){
        res[i] = a[i] - b[i];
    }
    return res;
}
std::vector<double> operator-(const std::vector<double> &a, double b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] + b;
    }
    return res;
}

std::vector<double> operator-=(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i){
        res[i] = a[i] - b[i];
    }
    return res;
}

double operator*(const std::vector<double> &a, const std::vector<double> &b) {
    double res=0;
    for (size_t i = 0, end = size(a); i < end; ++i){
        res += a[i] * b[i];
    }
    return res;
}

std::vector<double> operator*(double x, const std::vector<double> &a) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] * x;
    }
    return res;
}

std::vector<double> operator*(const std::vector<double> &a, double x) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] * x;
    }
    return res;
}

std::vector<double> operator/(const std::vector<double> &a, double x) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] / x;
    }
    return res;
}
std::vector<double> operator/(const std::vector<double> &a, const std::vector<double> &b) {
    const int n = size(a);
    std::vector<double> res(n);
    if (n != size(b)) {
        std::cout << "Dimensions of vectors must be equal! ('vec1/vec2' - element-wise division)" << std::endl;
        return res;
    }
    for (size_t i = 0; i < n; ++i) {
        res[i] = a[i]/b[i];
    }
    return res;
}

std::vector<double> operator/=(const std::vector<double> &x, const double a) {
    std::vector<double> res(size(x));
    for (size_t i = 0, end = size(x); i < end; ++i) {
        res[i] = x[i]/a;
    }
    return res;
}

bool operator==(const std::vector<double> &a, std::vector<double> &b) {
    if (size(a) != size(b)) {
        return false;
    }
    for (size_t i = 0, end=size(a); i < end; ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

// norms: first, second

double first_norm(const std::vector<double> &x) {
    double Norm = 0;
    for (size_t i = 0, end = size(x); i < end; ++i) {
        Norm += std::abs(x[i]);
    }
    return Norm;
}

double second_norm(const std::vector<double> &x) {
    return sqrt(x * x);
}

double energetic_norm(const std::vector<double> &x, const CSR &A) {
    return sqrt(x * (A * x));
}

//specific

std::vector<double> elWise_Mult(const std::vector<double> &a, const std::vector<double> &b) { //vect1 * vect2 == vect3, where vect3(i) == vect1(i) * vect2(i)
    const int n = size(a);
    std::vector<double> res(n);
    for (int i = 0; i < n; ++i) { 
        res[i] = a[i] * b[i];
    }
    return res;
}