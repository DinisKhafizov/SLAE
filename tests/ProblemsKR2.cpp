#include "PROBLEMS/MethodForKR2.hpp"
#include <fstream>

int main() {
    //first
    const int L = 17, a = 13, b2 = 62;
    const int N = L * L;
    const double tolerance = 0.000000000001;
    std::vector<double> vals(L * (5 * L - 2) - 2), b(N, 2), x_0(N);
    std::vector<int> rows(N + 1), cols(L * (5 * L - 2) - 2);
    vals[0] = b2;   cols[0] = 0;
    vals[1] = a;    cols[1] = 1;
    vals[2] = a;    cols[2] = L;
    rows[1] = 3;
    for (int i = 1 ; i < L; ++i) {
        vals[i*4 - 1] = a;  cols[i*4 - 1] = i - 1;
        vals[i*4] = b2;     cols[i*4] = i;
        vals[i*4 + 1] = a;  cols[i*4 + 1] = i + 1;
        vals[i*4 + 2] = a;  cols[i*4 + 2] = i + L;
        rows[i + 1] = rows[i] + 4;
    }
    for (int i = L; i < N - L; ++i) {
        vals[i*5 - L - 1] = a;  cols[i*5 - L - 1] = i - L;
        vals[i*5 - L] = a;      cols[i*5 - L] =  i - 1;
        vals[i*5 - L + 1] = b2; cols[i*5 - L + 1] = i;
        vals[i*5 - L + 2] = a;  cols[i*5 - L + 2] = i + 1;
        vals[i*5 - L + 3] = a;  cols[i*5 - L + 3] = i + L;
        rows[i+1] = rows[i] + 5;
    }
    for (int i = N - L; i < N - 1; ++i) {
        vals[i*4 + N - 2*L - 1] = a;    cols[i*4 + N - 2*L - 1] = i - L;
        vals[i*4 + N - 2*L] = a;    cols[i*4 + N - 2*L] = i - 1;
        vals[i*4 + N - 2*L + 1] = b2;       cols[i*4 + N - 2*L + 1] = i;
        vals[i*4 + N - 2*L + 2] = a;    cols[i*4 + N - 2*L + 2] = i + 1;
        rows[i+1] = rows[i] + 4;
    }
    vals[L * (5 * L - 2) - 5] = a; cols[L * (5 * L - 2) - 5] = N - 1 - L;
    vals[L * (5 * L - 2) - 4] = a; cols[L * (5 * L - 2) - 4] = N - 2;
    vals[L * (5 * L - 2) - 3] = b2; cols[L * (5 * L - 2) - 3] = N - 1;
    rows.back() = rows[N - 1] + 3;

    CSR A(vals, cols, rows, N);
    std::vector<double> res_si, res_si_opt, res_cheb, res_sym, res_sym_cheb, res1;
    std::ofstream fout1;
    fout1.open("Task11.txt");
    double counter = 0;
    std::vector<double> res = {counter, log(first_norm(A * x_0 - b))};  
    std::pair<double,double> lam = lambda(a, b2/2, L);
    res_si = SIM(A, x_0, b, tolerance, 1/lam.second); Write(res_si, "SIM.txt");
    res_si_opt = SIM(A, x_0, b, tolerance, 2/(lam.first + lam.second)); Write(res_si_opt, "SIM_W_OPT.txt");
    res_cheb = SIMwCA(A, x_0, b, tolerance, lam.first, lam.second, 3, res); Write(res_cheb, "SIM_W_CHEB.txt");
    res_sym = SGSM(A, x_0, b, tolerance); Write(res_sym, "SGS.txt");
    res_sym_cheb = SGSMwCA(A, x_0, b, tolerance, 0.5); Write(res_sym_cheb, "SGS_W_CHEB.txt");
    for (int i = 0; i < 3000; ++i) {
        res1 = SIMwCA_next(A, x_0, b, tolerance, lam.first, lam.second + i, 3);
        fout1 << res1[0] << ";" << res1[1] << "\n";
    }
    
    //second

    std::vector<double> vals2 = {14, 16, 18, 21}, vec_b(4, 4);
    std::vector<int> cols2 = {0, 1, 2, 3}, rows2 = {0, 1, 2, 3, 4};
    CSR A2(vals2, cols2, rows2, 4);
    double lambda_max = 21, lambda_min = 14;
    std::vecotr
}
