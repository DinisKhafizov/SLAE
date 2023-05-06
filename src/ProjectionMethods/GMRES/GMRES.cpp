#include "GMRES.hpp"


Hessenberg::Hessenberg(const CSR &A, const std::vector<double> &r_0, const int i_end ) {
    std::vector<double> v = r_0/second_norm(r_0), h(i_end + 1);
    const int N = A.GetN();
    double n;
    V.SetStringOnEnd(v);
    
    for (size_t j = 0; j < i_end; ++j) {    
        v = A * v;
        for (size_t k = 0, end = j + 1; k < end; ++k) {
            h[k] = v * V.getRow(k);
            v = v - h[k] * V.getRow(k);
        }
        n = second_norm(v);
        h[j + 1] = n;
        H.SetStringOnEnd(h);
        v = v/n;
        V.SetStringOnEnd(v);
    }
    H.transpose();
    V.transpose();
}

//i - how much extra iters you need
void Hessenberg::newIter(const CSR &A, const int i) { //If you getting another iterations for Hessenberg, you must give matrix as in constructor!
    const int N_now = V.GetN(), N = A.GetN();
    const int i_end = N_now + i;
    double n;
    std::vector<double> v = V.getCol(N_now - 1);
    std::vector<double> nulls(H.GetN());
    for (int j = 0; j < i; ++j) {
        H.SetStringOnEnd(nulls);
    }
    H.transpose();
    V.transpose();
    std::vector<double> h(H.GetN());
    if (size(h) == 1) {
        h.resize(2);
        H.GetN() = 2;
    }
    for (size_t j = N_now; j < i_end; ++j) {
        v = A * v;
        for (size_t k = 0, end = j; k < end; ++k) {
            h[k] = v * V.getRow(k);
            v = v - h[k] * V.getRow(k);
        }
        n = second_norm(v);
        h[j] = n;
        H.SetStringOnEnd(h);
        v = v/n;
        V.SetStringOnEnd(v);
    }
    H.transpose();
    V.transpose();
}

std::vector<double> Hessenberg::givens() { //don't need it actually
    const int N = H.GetN();
    double cos_t, sin_t, h1, h2, denom;
    std::vector<double> q(N + 1, 1);
    for (size_t i = 0; i < N; ++i) {
        h1 = H(i, i);
        h2 = H(i + 1, i);
        denom = sqrt(h1 * h1 + h2*h2);
        cos_t = h1/denom;
        sin_t = -h2/denom;
        q[i + 1] = q[i] * sin_t;
        q[i] *= cos_t;
        for (size_t j = i; j < N; ++j) {
            h1 = H(i, j);
            h2 = H(i + 1, j);
            H(i, j) = h1 * cos_t - h2 * sin_t;
            H(i + 1, j) = h1 * sin_t + h2 * cos_t;
        }
        H(i + 1, i) = 0;
    }
    return q;
}

std::vector<double> Hessenberg::givens_last_iter(std::vector<double> q) { //q - the first Q column except of first iteration
    const int N = H.GetN(), M = H.GetM(), r = size(rot), q_size = size(q);
    const int NmOne = N - 1;
    int k;
    double cos_t, sin_t, h1, h2, denom;
    if (r != 0) {
        rot.resize(r + 2);
        q.resize(q_size + 1);
        for (size_t i = 0, end = H.GetM() - 2; i < end; ++i) {
            k = i*2;
            h1 = H(i, NmOne);
            h2 = H(i + 1, NmOne);
            H(i, NmOne) = rot[k] * h1 - rot[k + 1] * h2;
            H(i + 1, NmOne) = rot[k + 1] * h1 + rot[k] * h2;
        }
        h1 = H(M - 2, NmOne);
        h2 = H(M - 1, NmOne);
        denom = sqrt(h1*h1 + h2*h2); 
        cos_t = h1/denom; 
        sin_t = -h2/denom;
        rot[r] = cos_t;
        rot[r + 1] = sin_t;
        q[q_size] = q[q_size - 1] * sin_t;
        q[q_size - 1] *= cos_t;
        H(M - 2, NmOne) = cos_t * h1 - sin_t * h2;
        H(M - 1, NmOne) = 0; //cos_t * h2 + sin_t * h1
    }
    else {
        rot.resize(2);
        h1 = H(0, 0);
        h2 = H(1, 0);
        denom = sqrt(h1*h1 + h2*h2);
        cos_t = h1/denom;
        sin_t = -h2/denom;
        rot[0] = cos_t;
        rot[1] = sin_t;
        q[0] = cos_t;
        q[1] = sin_t;
        H(0, 0) = h1 * cos_t - h2 * sin_t;
        H(1, 0) = 0; //h1 * sin_t + h2 * cos_t (rotation specifics)
    }
    return q;
}

Matrix Hessenberg::get_H() {
    return H;
}
Matrix Hessenberg::get_V() {
    return V;
}
Matrix Hessenberg::get_V_exc_lastcol() {
    Matrix B = V;
    B.transpose();
    B.deleteStringOnEnd();
    B.transpose();
    return B;
}


std::vector<double> GMRES(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, const double tolerance, int iters) {

    std::vector<double> r_0 = A * x_0 - b, x(size(x_0)), q, y;
    Hessenberg heis(A, r_0, iters);
    std::vector<double> q_all = heis.givens();
    double q_last = 1, norm_r0 = second_norm(r_0);
    for (size_t i = 0; i < iters; ++i) {
        if (std::abs(q_last * norm_r0) < tolerance) {
            return x;
        }
        q_last = q_all[i + 1];
        q.resize(i + 1);
        q[i] = q_all[i];
        x = x_0 - heis.get_V().partly_dot(GaussReverse(heis.get_H(), q * norm_r0, i + 1));
    }
    
    return x;
}

std::vector<double> GMRES(const CSR &A, const std::vector<double> &b, const std::vector<double> &x_0, const double tolerance) {
    
    const int N = A.GetN();
    std::vector<double> r_0 = A * x_0 - b, x(size(x_0)), q, y;
    Hessenberg heis(A, r_0, N);
    std::vector<double> q_all = heis.givens();
    double q_last = 1, norm_r0 = second_norm(r_0);
    for (size_t i = 0; i < N; ++i) {
        if (std::abs(q_last * norm_r0) < tolerance) {
            return x;
        }
        q_last = q_all[i + 1];
        q.resize(i + 1);
        q[i] = q_all[i];
        //y = GaussReverse(heis.get_H(), q, i + 1);
        x = x_0 - heis.get_V().partly_dot(GaussReverse(heis.get_H(), q * norm_r0, i + 1));
    }
    if (std::abs(q_last * norm_r0) < tolerance) {
        return x;
    }
        return GMRES(A, b, x, tolerance);
}