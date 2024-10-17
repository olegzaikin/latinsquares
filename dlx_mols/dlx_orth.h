#ifndef DLX_ORTH_H
#define DLX_ORTH_H

#include <vector>
#include <cassert>

#define latinsquare_t vector<vector<int>>

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

namespace DLX_orth {
	void cover(DLX_column *&c);
	void uncover(DLX_column *&c);
	void choose_c(DLX_column &h, DLX_column *&c);
	void square_to_DLX(DLX_column &root, const latinsquare_t square, vector<DLX_column*> &elements);
	void transversals_to_dlx(DLX_column &root, vector<vector<int>> &tvset, vector<DLX_column*> &elements);
	void find_all_transversals(int k, DLX_column &h, vector<DLX_column*> &ps, vector<vector<int>> &tvr);
	vector<vector<int>> find_tv_dlx(const latinsquare_t square);
	bool is_diag_latinsquare(const latinsquare_t square);
	bool is_latinsquare(const latinsquare_t square);
	vector<latinsquare_t> find_all_orth_mates(const latinsquare_t square);
};

#endif
