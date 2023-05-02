#ifndef HEY_ARNOLDI
#define HEY_ARNOLDI

#include "Matrixes/CSR/MatrixOnCSR.hpp"
#include "Matrixes/SMatrix/StandardMatrix.hpp"
#include "Vect/VectorOperations.hpp"

/*
std::pair<Matrix, Matrix> Arnoldi(const CSR &A, const std::vector<double> &r_0) {
    const int N = A.GetN();
    std::vector<double> V
}
*/

class Heisenberg {
    private:
        int i;
        Matrix H, V; //while filling H and V we get them transposed!
    public:
        Heisenberg(const CSR &A, const std::vector<double> &r_0, const int i_end) {
            std::vector<double> v = r_0/second_norm(r_0);
            const int N = A.GetN();
            double n;
            V.SetStringOnEnd(v);
            
            if (i_end >= N) {
                std::vector<double> h(N);
                for (size_t j = 0; j < N - 1; ++j) {    
                    v = A * v;
                    for (size_t k = 0, end = j + 1; k < end; ++k) {
                        h[k] = v * V.getRow(k);
                        v = v - h[k] * V.getRow(k);
                    }
                    n = second_norm(v);
                    h[j + 1] = n;
                    H.SetStringOnEnd(h);
                    v = v/n;
                    V.SetStringOnEnd(v);
                }
            }
            else {
                std::vector<double> h(i_end);
                for (size_t j = 0; j < i_end - 1; ++j) {    
                    v = A * v;
                    for (size_t k = 0, end = j + 1; k < end; ++k) {
                        h[k] = v * V.getRow(k);
                        v = v - h[k] * V.getRow(k);
                    }
                    n = second_norm(v);
                    h[j + 1] = n;
                    H.SetStringOnEnd(h);
                    v = v/n;
                    V.SetStringOnEnd(v);
                }
            }

            H.transpose();
            V.transpose();
        }
        /*
        void newIter(const int i) {
            std::vector<double> v = V.getCol(V.get)
            for (size_t j = 0; j < i; ++j) {
                v = A * v;
                for (size_t k = 0; k < j; ++k) {
                    h[k] = v * V.getRow(k);
                    v -= h[k] * V.getRow(k);
                }
                h[j] = second_norm(v);
                H.SetStringOnEnd(h);
                v /= second_norm(v);
            }
        }
        */

        Matrix get_H() {
            return H;
        }
        Matrix get_V() {
            return V;
        }

};

#endif