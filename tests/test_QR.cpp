#include <gtest/gtest.h>
#include "QR/Householder.hpp"

TEST(qr, householder1) {
    std::vector<double> a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    Matrix A(a, 4, 3);
    std::pair<Matrix, Matrix> QR = QR_Decomposition(A);
}
TEST(qr, householder2) {
    std::vector<double> a = {5, 0, 10, 13, 17, 2, 213, 32, 1235, 9942, 213, 1, 95, 21, 123, 4, 5, 1, 4, 41, 12, 123, 214, 21494, 1248, 124, 13, 9513};
    Matrix A(a, 7, 4);
    std::pair<Matrix, Matrix> QR = QR_Decomposition(A);
}