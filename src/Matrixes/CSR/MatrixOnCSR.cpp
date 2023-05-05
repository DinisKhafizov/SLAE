#include "MatrixOnCSR.hpp"

double CSR::operator()(int i, int j) const { 
	for(int k = rows[i]; k < rows[i + 1]; ++k) {
		if (cols[k] == j) {
			return values[k];
		}
	}
	return 0;
}

std::vector<double> CSR::operator*(const std::vector<double> &vec) const {
	std::vector<double> res(ROWS);	
	for (int i = 0; i < ROWS; ++i) {
		for (int j = rows[i]; j < rows[i + 1]; ++j) {
			res[i] += values[j] * vec[cols[j]];
		}
	}
	return res;
}
void CSR::operator *(double x) {
	for (int i = 0; i < size(values); ++i) {
		values[i] *= x;
	}
}

std::vector<double> CSR::getDiag() const { //only for square matrix!
	std::vector<double> res(COLS);
	for (size_t i = 0; i < COLS; ++i) {
		for (size_t k = rows[i]; k < rows[i + 1]; ++k) {
			if (cols[k] == i) { 
				res[i] = values[k];
			}
		}
	}
	return res;
}

int CSR::GetN() const {
	return ROWS;
}
std::vector<double> &CSR::GetVal() {
	return values;
}
std::vector<int> &CSR::GetRow()  {
	return rows;
}
std::vector<int> &CSR::GetCol() {
	return cols;
}
std::vector<double> CSR::GetVal() const {
	return values;
}
std::vector<int> CSR::GetRow()  const{
	return rows;
}
std::vector<int> CSR::GetCol()  const{
	return cols;
}
	