#include "PROBLEMS/MethodForProblem1.hpp"
#include <fstream>

int main() {
    const int L = 17, a = 13, b2 = 62;
    const int N = L * L;
    std::vector<double> rows(N + 1), vals(L * (5 * L - 2) - 2), cols(L * (5 * L - 2) - 2), b(N, 2);
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
    for (int i = L; i < N - 1 - L; ++i) {
        vals[i*5 - L - 1] = a;  cols[i*5 - L - 1] = i - L;
        vals[i*5 - L] = a;      cols[i*5 - L] =  i - 1;
        vals[i*5 - L + 1] = b2; cols[i*5 - L + 1] = i;
        vals[i*5 - L + 2] = a;  cols[i*5 - L + 2] = i + 1;
        vals[i*5 - L + 3] = a;  cols[i*5 - L + 3] = i + L;
        rows[i+1] = rows[i] + 5;
    }
    for (int i = N - 1 - L; i < N - 1; ++i) {
        vals[i*4 + N - 2*L - 2] = a;    cols[i*4 + N - 2*L - 2] = i - L;
        vals[i*4 + N - 2*L - 1] = a;    cols[i*4 + N - 2*L - 1] = i - 1;
        vals[i*4 + N - 2*L] = b2;       cols[i*4 + N - 2*L] = i;
        vals[i*4 + N - 2*L + 1] = a;    cols[i*4 + N - 2*L + 1] = i + 1;
        rows[i+1] = rows[i] + 4;
    }
    vals[]
}
