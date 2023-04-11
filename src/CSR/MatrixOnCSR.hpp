#ifndef CSRMATRIX
#define CSRMATRIX

#include<iostream>
#include<vector>

class CSR {
	int ROWS;
	int COLS;
	std::vector<double> values;
	std::vector<int> cols;
	std::vector<int> rows;
	std::vector<double> dopEls;
public:

	CSR(const int n, const int m): ROWS{m}, COLS{n}, values(m * n), cols(0), rows(1, 0) 
	{} 

	CSR(const std::vector<double> &val, const std::vector<int> &col, const std::vector<int> &row, const int Ncols)
	:values{val}, cols{col}, rows{row}, ROWS{size(row) - 1}, COLS{Ncols} 
	{}

	double operator()(int i, int j) const { 
		for(int k = rows[i]; k < rows[i + 1]; ++k) {
			if (cols[k] == j) {
				return values[k];
			}
		}
		return 0;
	}

	std::vector<double> operator*(const std::vector<double> &vec) const {
	    std::vector<double> res(ROWS);	
		for (int i = 0; i < ROWS; ++i) {
			for (int j = rows[i]; j < rows[i + 1]; ++j) {
				res[i] += values[j] * vec[cols[j]];
			}
		}
		return res;
	}
	void operator *(double x) {
		for (int i = 0; i < size(values); ++i) {
			values[i] *= x;
		}
	}
	int GetN() const {
		return ROWS;
	}
	std::vector<double> &GetVal() {
		return values;
	}
	std::vector<int> &GetRow() {
		return rows;
	}
	std::vector<int> &GetCol() {
		return cols;
	}
	
};

#endif //CSRMATRIX