#ifndef GMRES_meth
#define GMRES_meth

#include "Matrixes/CSR/MatrixOnCSR.hpp"
#include "Matrixes/SMatrix/StandardMatrix.hpp"
#include "Vect/VectorOperations.hpp"
#include "Gauss/GaussReverse.hpp"

class Hessenberg {
    private:
        Matrix H, V; //while filling H and V we get them transposed!
        std::vector<double> rot; // begins with cos
    public:
        Hessenberg(const CSR &A, const std::vector<double> &r_0, const int i_end =1);
        //i - how much extra iters you need
        void newIter(const CSR &A, const int i = 1);

        std::vector<double> givens();

        std::vector<double> givens_last_iter(std::vector<double> q = {1, 1});

        Matrix get_H();
        Matrix get_V();
        Matrix get_V_exc_lastcol();

};


std::vector<double> GMRES(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, const double tolerance, int iters);

std::vector<double> GMRES(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, const double tolerance);

#endif