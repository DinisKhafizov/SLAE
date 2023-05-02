//#ifndef QR_MATRIX
//#define QR_MATRIX

#include <iostream>
#include <vector>
//#include "SMatrix/StandardMatrix.hpp"

//::pair<Matrix, Matrix> 

class Matrix {
private:
    std::vector<double> A;
    int M, N;
public:
    Matrix(const int m, const int n) {
        std::vector<double> B(m * n);
        A = B;
        M = m;
        N = n;
    }
    Matrix(const std::vector<double> &a, const int m, const int n) { //N - кол-во столбцов
        A = a;
        N = n;
        M = m;
    }
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







Matrix QR_Decomposition(Matrix &A) {

    std::vector<double> x, v, P_first_str, x_other, vall;

    double v_quad = 0, v_x = 0, mod_x = 0, r_ii = 0;

    int size_v, size_vall;

    int N = A.GetN(), M = A.GetM();

    Matrix R = A, Q(M, M), v_all(0, N); 

    //                                                      COUNTING R

    for (int i = 0; i < N; ++i) {
        x = R.GetCol_WOut_Els(i, i); //column "x" from matrix (P_i * ... * P_1) * A without first "i" strings

        for (int k = 0; k < size(x); ++k) {
            mod_x += x[k] * x[k];    //get module of x to make vector v
        }
        
        v = x;                       // creating v with adding module x at first (of v) coordinate
        if (x[0] > 0) {
            v[0] += mod_x;
            }
        else {
            v[0] -= mod_x;
        }

        size_v = size(v);
        size_vall = size(vall);
        vall = {0};
        vall.resize(N);

        for (int k = 0; k < size_v; ++k) {
            vall[k] = v[k];
        }                            // adding vector v to matrix v_all that saves v at each iteration

        for (int k = size(v); k < N; ++k) {
            vall.push_back(0);
        }

        v_all.SetStringOnEnd(vall);

        for (int k = 0; k < size(v); ++k) {
            v_quad += v[k] * v[k];   //counting scalar (v * v)
        }

        r_ii += (1 - v[0] * v[0] / v_quad) * x[0];
        for (int k = 1; k < size(v); ++k) {
            r_ii += (2 * v[0] * v[k] / v_quad) * x[k];
        }                             // counting the el that would be on diagonal on (i, i) place

        R(i, i) = r_ii;               // giving the R this el

        for (int k = i + 1; k < M; ++k) {
              R(k, i) = 0;            // with following nulling of els which are lower than diag el
        }
        for (int j = i + 1; j < N; ++j) { //here is counting els after i column
            
            x_other = R.GetCol_WOut_Els(j, i); 

            for (int k = 0; k < size(x_other); ++k) { //counting scalar x_other on v
                v_x += v[k] * x_other[k];
            }

            for (int k = i; k < size(x_other); ++k) { //writing new els in our R
                R(k, j) = x_other[k] - 2 * v_x / v_quad * v[0];
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        std::cout << std::endl;
        for (int j = 0; j < 3; ++j) {
            std::cout << v_all(i, j) << " ";
        }
    }
    std::cout << std::endl;
    return R;

    //                                                      COUNTING Q

    


} 


int main() {
    std::vector<double> a = {1, 2, 3, 4, 5, 6, 7, 8, 10};
    int n = 3, m = 3;
    Matrix A(a, m, n), res(m, n);
    res = QR_Decomposition(A);
    std::cout << "===============" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << std::endl;
        for (int j = 0; j < 3; ++j) {
            std::cout << res(i, j) << " ";
        }
    }

    return 0;
}

//#endif 
