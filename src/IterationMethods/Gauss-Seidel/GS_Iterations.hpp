#ifndef ITERS
#define ITERS

#include<iostream>
#include"Vect/VectorOperations.hpp"
#include"Matrixes/CSR/MatrixOnCSR.hpp"
#include "Vect/VectorOperations.hpp"

//Комменты:

//11.04: улучшил итерации сверху вниз и наоборот. Лишнюю итерацию с прибавлением x[i] += x[i] * diag[i] убрал
//посредством зануления диагональных элементов у матрицы А (у вектора values CSR матрицы). Других идей по оптимизации итераций пока нет.

//05.05: матрица теперь передается по ссылке, const. Более не зануляю диагональные элементы, производится проверка на диагональный элемент,
//пока идет итерация по строке A. Оптимизировано.

std::vector<double> TopDownIteration(const CSR &A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x);

std::vector<double> DownUpIteration(const CSR &A, const std::vector<double> &diag, const std::vector<double> &b, std::vector<double> x);

#endif