#include <gtest/gtest.h>
#include "SolverTridiag/Solver.hpp"

TEST(solver, slae1) {
    std::vector<double> a = {1., 1.};
    std::vector<double> b = {4., 4., 4.};
    std::vector<double> c = {1., 1.};
    std::vector<double> f = {28., 7., 56.};

    TridiagonalMatrix A(a, b, c);
    std::vector<double> x;
    x = Progonka(A, f);

    ASSERT_DOUBLE_EQ(x[0], 8.);
    ASSERT_DOUBLE_EQ(x[1], -4.);
    ASSERT_DOUBLE_EQ(x[2], 15.);
}

TEST(solver, slae2) {
    std::vector<double> a = {5., 1.};
    std::vector<double> b = {2., 4., -3.};
    std::vector<double> c = {-1., 2.};
    std::vector<double> f = {3., 6., 2.};

    TridiagonalMatrix A(a, b, c);
    std::vector<double> x;
    x = Progonka(A, f);

    ASSERT_DOUBLE_EQ(x[0], 64./43.) ;
    ASSERT_DOUBLE_EQ(x[1], -1./43.);
    ASSERT_DOUBLE_EQ(x[2], -29./43.);
}

TEST(solver, slae3) {
    std::vector<double> a = {1., 3., 6., 2., 3., 11., 10., 1., 1.};
    std::vector<double> b = {3.3, 5.2, 6.37, 12.03, 3., 5.6, 95., 102., 2., 1.};
    std::vector<double> c = {2., 1., 2., 5., 1., 1., 50., 30., 0.5};
    std::vector<double> f = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    TridiagonalMatrix A(a, b, c);
    std::vector<double> x;
    x = Progonka(A, f);
    std::vector<double> x_ref = {362572824060./1551876085453., 355385766055./3103752170906.,
                                    265300269650./1551876085453., -335582640650./1551876085453.,
                                    1598853453829./3103752170906., -350477627981./3103752170906.,
                                    1349332630563./15518760854530., -186783139020./1551876085453.,
                                    1925442363493./4655628256359., 2730185892866./4655628256359.};
    for (std::size_t i = 0; i < 10; ++i) {
        ASSERT_DOUBLE_EQ(x[i], x_ref[i]);
    }
}
