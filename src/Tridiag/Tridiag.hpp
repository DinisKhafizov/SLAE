#ifndef THREE_DIAG_SOLVER_TRIDIAG_SOLVER_H
#define THREE_DIAG_SOLVER_TRIDIAG_SOLVER_H


#include <iostream>
#include <vector>
#include <string>
#include <fstream>


class TridiagonalMatrix 
{
private:
    int n;
    std::vector<double> A, B, C;
public:

//Пока еще разбираюсь с файловой системой
/*
    TridiagonalMatrix(std::string path){
        setlocale(LC_ALL, "ru");
        std::ifstream fin;
        fin.open(path);
        if (!fin.is_open()) {std::cout << "Ошибка чтения файла!" << std::endl;}
        fin >> n; 
        A.reserve(n - 1);
        B.reserve(n);
        C.reserve(n - 1);
        double flag;
        for (int i = 0; i < 3 * n - 2; ++i) {
            fin >> flag;
            if (i < n - 1) { C.push_back(flag); }
            else if (i < 2 * n - 1) { B.push_back(flag); }
            else if (i < 3 * n - 2) { A.push_back(flag); }
        }
        fin.close();
    
    }
*/
    TridiagonalMatrix(const std::vector<double> &A1, const std::vector<double> &B1,const std::vector<double> &C1){
        A = A1;
        B = B1;
        C = C1;
        n = size(B);
    }

    int GetN() {
        return n;
    }
    double GetA(int i) {
        return A[i];
    }
    double GetB(int i) {
        return B[i];
    }
    double GetC(int i) {
        return C[i];
    }
};

#endif //THREE_DIAG_SOLVER_TRIDIAG_SOLVER_H