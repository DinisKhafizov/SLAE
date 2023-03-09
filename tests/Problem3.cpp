#include "Work/MethodForProblem3.hpp"

void Problem_3() { 
	std::ofstream fout;
	fout.open("Problem3.txt");

	std::vector<double> a = {10, -2, 6, 3, 8, -1, 1, 2, 1}, b = {1, 2, 3}, x_0 = {0, 0, 0};
	std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
	double tolerance = 0.000000000001, tau = 0.001;
	int res;
	int Ncols = 3;
	CSR A(a, cols, rows, Ncols);
	for (int i = 1; i < 100; ++i) {
		res = SimpleIterationMethod(A, x_0, b, tolerance, tau * double(i));
		fout << res << ";" << tau * i << "\n";
	}
	fout.close();
}

int main() {
    Problem_3();
    return 0;
}