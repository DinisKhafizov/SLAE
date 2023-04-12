#ifndef VECTOR_OP
#define VECTOR_OP

#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (int i =0; i < size(a); ++i){
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (int i =0; i < size(a); ++i){
        res[i] = a[i] - b[i];
    }
    return res;
}


/*
Matrix operator*(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(size(a) * size(b));
    int k, sB = size(b), sA = size(a);
    for (int i = 0; i < size(a); ++i){
        k = i * i * sB;
        for (int j = 0; i < size(b); ++j){
            res[k + j] = a[k] * b[j];
        }
    }
    Matrix A(res, sA, sB);
    return A;
}
*/

std::vector<double> elWiseMult(const std::vector<double> &a, const std::vector<double> &b) {
    const int n = size(a);
    std::vector<double> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] * b[i];
    }
    return res;
}

double operator*(const std::vector<double> &a, const std::vector<double> &b) {
    double res=0;
    for (int i = 0; i < size(a); ++i){
        res += a[i] * b[i];
    }
    return res;
}

std::vector<double> operator*(double x, const std::vector<double> &a) {
    std::vector<double> res(size(a));
    for (int i = 0; i < size(a); ++i) {
        res[i] = a[i] * x;
    }
    return res;
}

std::vector<double> operator*(const std::vector<double> &a, double x) {
    std::vector<double> res(size(a));
    for (int i = 0; i < size(a); ++i) {
        res[i] = a[i] * x;
    }
    return res;
}

std::vector<double> operator/(const std::vector<double> &a, double x) {
    std::vector<double> res(size(a));
    for (int i = 0; i < size(a); ++i) {
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
    for (int i = 0; i < n; ++i) {
        res[i] = a[i]/b[i];
    }
    return res;
}

double first_norm(const std::vector<double> &x) {
    double Norm = 0;
    for (int i = 0; i < size(x); ++i) {
        Norm += std::abs(x[i]);
    }
    return Norm;
}

bool vect_equality(const std::vector<double> &a, std::vector<double> &b) {
    if (size(a) != size(b)) {
        return false;
    }
    for (int i = 0; i < size(a); ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

void head(std::vector<double> x, const int num = 10, int cols = 0) {
    for (int i = 0; i < num; ++i) {
        std::cout << "| ";
        for (int j = 0; j < cols; ++j){
            std::cout << x[i*cols + j] << " | ";
        }
        std::cout << std::endl << "--------" << std::endl;
    }
}


/*
std::vector<double> addition(const std::vector<double> &a,const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (int i = 0; i < size(a); ++i){
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> difference(const std::vector<double> &a,const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (int i = 0; i < size(a); ++i){
        res[i] = a[i] - b[i];
    }
    return res;
}

double mult(const std::vector<double> &vect_str, const std::vector<double> &vect_col) {
    double res = 0;
    for (int i = 0; i < size(vect_col); ++i) {
        res += vect_str[i] * vect_col[i];
    }
    return res;
}

std::vector<double> mult(const std::vector<double> &vect, double x) {
    std::vector<double> res(size(vect));
    for (int i = 0; i < size(vect); ++i) {
        res[i] = vect[i] * x;
    }
    return res;
}
std::vector<double> mult(double x, const std::vector<double> &vect) {
    std::vector<double> res(size(vect));
    for (int i = 0; i < size(vect); ++i) {
        res[i] = vect[i] * x;
    }
    return res;
}

std::vector<double> div(const std::vector<double> &vect, double x) {
    std::vector<double> res(size(vect));
    for (int i = 0; i < size(vect); ++i) {
        res[i] = vect[i] / x;
    }
    return res;
}


*/
#endif