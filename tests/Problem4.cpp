#include "Work/MethodsForProblem4.hpp"


void Problem_4() {
	std::ofstream fout_Jac, fout_Seid, fout_Simp;
	fout_Jac.open("Problem4_Jacob.txt");
	fout_Seid.open("Problem4_Seid.txt");
	fout_Simp.open("Problem4_Simp.txt");

	std::vector<double> a = {12, 7, 112, 7, 15367, 401, 112, 401, 7182}, b = {75, 438, 139}, x_0={0, 0, 7};
	std::vector<int> rows = {0, 3, 6, 9}, cols = {0, 1, 2, 0, 1, 2, 0, 1, 2};
	double tolerance = 0.00000001, tau = 0.00005;
	std::pair<std::vector<double>, std::vector<int>> res_Jac(std::vector<double>(0), std::vector<int>(0)), 
	res_Seid(std::vector<double>(0), std::vector<int>(0)), res_Simp(std::vector<double>(0), std::vector<int>(0));
	int Ncols = 3;
	CSR A(a, cols, rows, Ncols);
	res_Jac = JacobiMethod(A, x_0, b, tolerance);
	res_Seid = GaussSeidelMethod(A, x_0, b, tolerance);
	res_Simp = SimpleIterationMethod(A, x_0, b, tolerance, tau);
	for (int i = 0; i < size(res_Jac.first); ++i) {
		fout_Jac << res_Jac.first[i] << ";" << res_Jac.second[i] << std::endl;
	}
	for (int i = 0; i < size(res_Seid.first); ++i) {
		fout_Seid << res_Seid.first[i] << ";" << res_Seid.second[i] << std::endl;
	}
	for (int i = 0; i < size(res_Simp.first); ++i) {
		fout_Simp << res_Simp.first[i] << ";" << res_Simp.second[i] << std::endl;
	}

	fout_Jac.close();
	fout_Seid.close();
	fout_Simp.close();
}


int main() {
	Problem_4();
	return 0;
}
