#include <gtest/gtest.h>
#include "IterationMethods/Gauss-Seidel/GS.hpp"
#include "IterationMethods/Gauss-Seidel/SGS.hpp"
#include "IterationMethods/Gauss-Seidel/SGSwCA.hpp"
#include "IterationMethods/Gradient/CG.hpp"
#include "IterationMethods/Gradient/SGD.hpp"
#include "IterationMethods/Jacobi/Jacobi.hpp"
#include "IterationMethods/SimpleIteration/SI.hpp"
#include "IterationMethods/SimpleIteration/SIwCA.hpp"
#include "IterationMethods/SOR/SOR.hpp"

TEST(gs, gs1) {
    std::vector<double> b = {3, 9, 20 ,14}, x_0={0, 0, 0, 0}, vals = {1, 2, 3, 5, 2, 1, 10, 13, 60, 1, 2, 17, 19}, res = {-12.279, 3.718, 1.564, 0.63}, res1;
    std::vector<int> rows = {0, 4, 6, 10, 13}, cols = {0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 2, 3};
    CSR A(vals, cols, rows, 4);
    res1 = GSM(A, x_0, b, 0.000001);
    for (int i = 0; i < 4; ++i) {
        ASSERT_NEAR(res[i], res1[i], 0.001);
    }

}
TEST(sgs, sgs1){ 
    std::vector<double> b = {12, 51, 3, 123}, x_0={0, 0, 0, 0}, vals = {5, -1, -2, -1, 4, 1, 1, 3, -1, -2, -1, 4}, res = {23.718, 15.987, 10.772, 45.302}, res1;
    std::vector<int> rows = {0, 3, 6, 9, 12}, cols = {0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3};
    CSR A(vals, cols, rows, 4);
    res1 = SGSM(A, x_0, b, 0.000001);
    for (int i =0; i < 4; ++i) {
        std::cout << res1[i] << std::endl;
    }
    for (int i = 0; i < 4; ++i) { 
        ASSERT_NEAR(res[i], res1[i], 0.9);
    }

}
TEST(sgswca, sgswca1) {
    std::vector<double> b = {12, 51, 3, 123}, x_0={0, 0, 0, 0}, vals = {5, -1, -2, -1, 4, 1, 1, 3, -1, -2, -1, 4}, res = {23.718, 15.987, 10.772, 45.302}, res1;
    std::vector<int> rows = {0, 3, 6, 9, 12}, cols = {0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3};
    double lambda_max = 1.1;
    CSR A(vals, cols, rows, 4);
    res1 = SGSMwCA(A, x_0, b, 0.000001, lambda_max);
    for(int i = 0; i < 4; ++i) {
        std::cout << res1[i] << std::endl;
    }
    for (int i = 0; i < 4; ++i) { 
        ASSERT_NEAR(res[i], res1[i], lambda_max);
    }

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
