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
    
    double operator ()(int i, int j) const {
        return A[i * N + j];
    }
    double &operator ()(int i, int j) {
        return A[i * N + j];
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
    std::vector<double> GetMinor(int k) const {
        std::vector<double> vec;
        vec.reserve((N - k) * (N - k));
        for (int i = 0; i < (N - i) * (N - i); ++i) {
           // vec[i] = A[]
        }
    }

    void SetStringOnEnd(const std::vector<double> &newstr) {
        A.resize(M * N + N);
        for (int i = 0; i < N; ++i) {
            A[M * N + i] = newstr[i];
        }
        M += 1;
    }

    int GetN() {
        return(N);
    }

    int GetM() {
        return(M);
    }
};

#endif