#ifndef SMATRIX
#define SMATRIX

#include <iostream>
#include <vector>

class Matrix {
private:
    std::vector<double> A;
    int M, N;
public:
    Matrix(const int m, const int n): A{std::vector<double>(m* n)}, M{m}, N{n} {}
    Matrix(const std::vector<double> &a, const int m, const int n): A{a}, N{n}, M{m} {}
    Matrix(): N{0}, M{0} {}

    double operator ()(int i, int j) const;
    double &operator ()(int i, int j);

    std::vector<double> operator*(const std::vector<double> &x) const;

    std::vector<double> GetCol_WOut_Els(int j, int el) const;

	std::vector<double> GetStr_WOut_Els(int i, int el) const;

    std::vector<double> getRow(const int i);

    std::vector<double> getCol(const int j);


    void transpose();


    Matrix transpose(const int what);

    void SetStringOnEnd(const std::vector<double> &newstr);

    void deleteStringOnEnd();

    int &GetN();
    int &GetM();

    std::vector<double> &get_vals();

    void show();
};

#endif