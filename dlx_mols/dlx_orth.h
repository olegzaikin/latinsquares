#ifndef DLX_ORTH_H
#define DLX_ORTH_H

#include <vector>
#include <string>
#include <set>
#include <cassert>

#define row_t vector<unsigned short int>
#define latinsquare_t vector<row_t>
#define transversal_t vector<int>

using namespace std;

struct DLX_column {
	int size;
	int column_number;
	int row_id;

	DLX_column *Left;
	DLX_column *Right;
	DLX_column *Up;
	DLX_column *Down;
	DLX_column *Column;
};

struct LS_result {
	vector<latinsquare_t> diag_orth_mates;
	unsigned transv;
	unsigned diag_transv;
};

namespace DLX_orth {
	void cover(DLX_column *&c);
	void uncover(DLX_column *&c);
	void choose_c(DLX_column &h, DLX_column *&c);
	void square_to_DLX(DLX_column &root, const latinsquare_t square, vector<DLX_column*> &elements, const bool is_diag);
	void transversals_to_dlx(DLX_column &root, vector<vector<int>> &tvset, vector<DLX_column*> &elements);
	void find_transversals(int k, DLX_column &h, vector<DLX_column*> &ps, vector<transversal_t> &tvr);
	vector<vector<int>> find_tv_dlx(const latinsquare_t square, const bool is_diag);
	bool is_diag_latinsquare(const latinsquare_t square);
	bool is_latinsquare(const latinsquare_t square);
	unsigned calc_orth_char(const latinsquare_t dls1, const latinsquare_t dls2);
	bool is_pseudotriple(const latinsquare_t square1, const latinsquare_t square2, 
		const latinsquare_t square3, const unsigned n, const unsigned orth_char);
	string square_to_str(const latinsquare_t square);
	vector<latinsquare_t> find_orth_mates(const latinsquare_t square);
	LS_result find_diag_orth_mates(const latinsquare_t square);
};

#endif
