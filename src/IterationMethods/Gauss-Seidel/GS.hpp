#ifndef GAUSS_SEID_METH
#define GAUSS_SEID_METH

#include <iostream>
#include <vector>
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"
#include "GS_Iterations.hpp"

//Комменты:

//05.05: матрица теперь передается по ссылке, const. Более не зануляю диагональные элементы, производится проверка на диагональный элемент,
//пока идет итерация по строке A. Оптимизировано.

std::vector<double> GSM(const CSR &A, const std::vector<double> &x_0, const std::vector<double> &b, const double tolerance);


#endif