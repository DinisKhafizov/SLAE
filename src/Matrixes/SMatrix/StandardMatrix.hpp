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

    double operator ()(int i, int j) const {
        return A[i * N + j];
    }
    double &operator ()(int i, int j) {
        return A[i * N + j];
    }

    std::vector<double> operator*(const std::vector<double> &x) const {
        std::vector<double> res(M);
        int k;
        for (size_t i = 0; i < M; ++i) {
            k = i * N;
            for (size_t j = 0; j < N; ++j) {
                res[i] += A[k + j] * x[j];
            }
        }
        return res;
    }

    std::vector<double> GetCol_WOut_Els(int j, int el) const	{
		std::vector<double> vec(M - el);
		for (int i = 0 + el; i < M; ++i){
			vec[i - el] = A[i * N + j];
		}
		return vec;
    }

	std::vector<double> GetStr_WOut_Els(int i, int el) const	{
		std::vector<double> vec;
		vec.reserve(N);
		for (int j = 0 + el; j < N; ++j){
			vec[j] = A[i * M + j];
		}
		return vec;
    }

    std::vector<double> getRow(const int i) {
        std::vector<double> vec(N);
        const int k = i * N;
        for (int j = 0; j < N; ++j) {
            vec[j] = A[k + j];
        }
        return vec;
    }

    std::vector<double> getCol(const int j) {
        std::vector<double> vec(M);
        for (int i = 0; i < M; ++i) {
            vec[i] = A[j - 1 + i*N];
        }
        return vec;
    }


    void transpose() {
        std::vector<double> res(M * N);
        int k;
        for (size_t i = 0; i < N; ++i) {
            k = i * M;
            for (size_t j = 0; j < M; ++j) {
                res[j + k] = A[i + j*N];
            }
        }
        A = res;
        k = N;
        N = M;
        M = k;
    }


    Matrix transpose(const int what){
        std::vector<double> res(M * N);
        int k;
        for (size_t i = 0; i < N; ++i) {
            k = i * M;
            for (size_t j = 0; j < M; ++j) {
                res[j + k] = A[i + j*N];
            }
        }
        Matrix B(res, N, M);
        return B;
    }

    void SetStringOnEnd(const std::vector<double> &newstr) {
        const int n = size(newstr), m = size(A);
        A.resize(m + n);
        for (int i = 0; i < n; ++i) {
            A[m + i] = newstr[i];
        }
        if (M == 0) {
            M = 1;
        }
        else {
            M += 1;
        }
        if (N == 0) {
            N = n;
        }
    }

    void deleteStringOnEnd() {
        for (size_t i = 0; i < N; ++i) {
            A.pop_back();
        }
        M -= 1;
    }

    int &GetN() {
        return(N);
    }

    int &GetM() {
        return(M);
    }

    std::vector<double> &get_vals() {
        return A;
    }

    void show() {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << A[i*N + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};

#endif