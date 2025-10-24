#include <algorithm>
#include <set>
#include <cassert>

#include "dlx_orth.h"

// Check whether a given square is a Latin square:
bool DLX_orth::is_latinsquare(const latinsquare_t square) {
	const unsigned n = square.size();
	assert(n > 0);
	for (unsigned i = 0; i < n; i++) {
		assert(square[i].size() == n);
		// Check the row's elements uniqueness:
		std::set<int> row;
		for (auto &x : square[i]) row.insert(x);
		if (row.size() < n) return false;
		// For each column, check whether elements are unique:
		std::set<int> column;
		for (unsigned j = 0; j < n; j++) {
			column.insert(square[j][i]);
		}
		if (column.size() < n) return false;
	}
	return true;
}

// Check whether a given square is a diagonal Latin square:
bool DLX_orth::is_diag_latinsquare(const latinsquare_t square) {
	const unsigned n = square.size();
	assert(n > 0);
	if (not is_latinsquare(square)) return false;
	std::set<int> main_diag;
	std::set<int> main_antidiag;
	for (unsigned i = 0; i < n; i++) {
		main_diag.insert(square[i][i]);
		main_antidiag.insert(square[i][n - i - 1]);
	}
	if (main_diag.size() < n or main_antidiag.size() < n) return false;
	return true;
}

void DLX_orth::cover(DLX_column *&c) {
	//cout << "Covered " << c->column_number << endl;
	c->Right->Left = c->Left;
	c->Left->Right = c->Right;

	DLX_column *i;
	DLX_column *j;
	i = c->Down;
	while (i != c) {
		j = i->Right;
		while (j != i) {
			j->Down->Up = j->Up;
			j->Up->Down = j->Down;
			//	cout << "covered element " << j->row_id << " in column " << j->column_number << endl;
			j->Column->size--;
			//if (j->Column->size < 0) {
			//	cout << "trouble" << endl;
			//}
			j = j->Right;
		}
		i = i->Down;
	}
}

void DLX_orth::uncover(DLX_column *&c) {
	//cout << "Uncovered " << c->column_number << endl;
	DLX_column *i;
	DLX_column *j;
	i = c->Up;
	while (i != c) {
		j = i->Left;
		while (j != i) {
			j->Column->size++;

			j->Down->Up = j;
			j->Up->Down = j;

			j = j->Left;
		}
		i = i->Up;
	}
	c->Right->Left = c;
	c->Left->Right = c;
}

void DLX_orth::choose_c(DLX_column &h, DLX_column *&c) {
	DLX_column * j;

	j = h.Right;
	int min = j->size;
	c = j;
	while (j != &h) {
		if (j->size < min) {
			c = j;
			min = j->size;
		}
		j = j->Right;
	}
}

void DLX_orth::square_to_DLX(DLX_column &root, const latinsquare_t square, vector<DLX_column*> &elements, const bool is_diag) {
	const unsigned n = square.size();
	assert(n > 0);
	root.Up = NULL;
	root.Down = NULL;
	root.Column = NULL;
	root.row_id = -1;
	root.size = -1;
	//	root.column_number= -1;
	elements.push_back(&root);
	vector<DLX_column *> columns;
	DLX_column * lastleft = &root;
	// first n - row number
	// n to 2n - column number
	//2n to 3n - value
	//3n+1 - diag
	//3n+2 - antidiag
	int dlx_columns_num;
	if (is_diag) dlx_columns_num = 3*n+2;
	else dlx_columns_num = 3*n;
	for (int i = 0; i < dlx_columns_num; i++) {
		DLX_column *ct;
		ct = new (DLX_column);
		//	ct->column_number = i;
		ct->Down = ct;
		ct->Up = ct;
		ct->size = 0;
		ct->row_id = 0;
		ct->Column = ct;
		ct->Left = lastleft;
		lastleft->Right = ct;
		lastleft = ct;
		columns.push_back(ct);
		elements.push_back(ct);
	}
	lastleft->Right = &root;
	root.Left = lastleft;
	for (int i = 0; i < n; i++) {
		assert(square[i].size() == n);
		for (int j = 0; j < n; j++) {
			vector<DLX_column *> tvrow;
			DLX_column *ctve;
			ctve = new (DLX_column);

			ctve->Column = columns[i];
			ctve->Column->size++;
			ctve->Down = columns[i];
			ctve->Up = columns[i]->Up;
			ctve->Up->Down = ctve;
			ctve->Down->Up = ctve;
			ctve->row_id = i*n +j;
			//	ctve->column_number = k;
			ctve->size = -10;
			elements.push_back(ctve);
			tvrow.push_back(ctve);

			ctve = new (DLX_column);
			//column corresponds to characteristic vector of LS or smth of that kind

			ctve->Column = columns[n + j];
			ctve->Column->size++;
			ctve->Down = columns[n + j];
			ctve->Up = columns[n + j]->Up;
			ctve->Up->Down = ctve;
			ctve->Down->Up = ctve;
			ctve->row_id = i*n + j;
			//	ctve->column_number = k;
			ctve->size = -10;
			elements.push_back(ctve);
			tvrow.push_back(ctve);

			ctve = new (DLX_column);
			//column corresponds to characteristic vector of LS or smth of that kind

			ctve->Column = columns[2*n + square[i][j]];
			ctve->Column->size++;
			ctve->Down = columns[2 * n + square[i][j]];
			ctve->Up = columns[2 * n + square[i][j]]->Up;
			ctve->Up->Down = ctve;
			ctve->Down->Up = ctve;
			ctve->row_id = i*n + j;
			//	ctve->column_number = k;
			ctve->size = -10;
			elements.push_back(ctve);
			tvrow.push_back(ctve);

			if (is_diag) {
				if (i == j) {
					ctve = new (DLX_column);
					ctve->Column = columns[3 * n ];
					ctve->Column->size++;
					ctve->Down = columns[3 * n ];
					ctve->Up = columns[3 * n]->Up;
					ctve->Up->Down = ctve;
					ctve->Down->Up = ctve;
					ctve->row_id = i*n + j;
					//	ctve->column_number = k;
					ctve->size = -10;
					elements.push_back(ctve);
					tvrow.push_back(ctve);
				}
				if (i == (n - j - 1)) {
					ctve = new (DLX_column);
					ctve->Column = columns[3 * n+1];
					ctve->Column->size++;
					ctve->Down = columns[3 * n+1];
					ctve->Up = columns[3 * n+1]->Up;
					ctve->Up->Down = ctve;
					ctve->Down->Up = ctve;
					ctve->row_id = i*n + j;
					//	ctve->column_number = k;
					ctve->size = -10;
					elements.push_back(ctve);
					tvrow.push_back(ctve);
				}	
			}	

			for (int j = 0; j < tvrow.size() - 1; j++) {
				tvrow[j]->Right = tvrow[j + 1];
				tvrow[j]->Right->Left = tvrow[j];
			}
			tvrow[tvrow.size() - 1]->Right = tvrow[0];
			tvrow[0]->Left = tvrow[tvrow.size() - 1];

		}
	}
	DLX_column *pr = &root;

}

void DLX_orth::transversals_to_dlx(DLX_column &root, vector<vector<int>> &tvset, vector<DLX_column*> &elements) {
	int n = tvset[0].size();
	root.Up = NULL;
	root.Down = NULL;
	root.Column = NULL;
	root.row_id = -1;
	root.size = -1;
	//	root.column_number= -1;
	elements.push_back(&root);
	vector<DLX_column *> columns;
	DLX_column * lastleft = &root;
	for (int i = 0; i < n* n; i++) {
		DLX_column *ct;
		ct = new (DLX_column);
		//	ct->column_number = i;
		ct->Down = ct;
		ct->Up = ct;
		ct->size = 0;
		ct->row_id = 0;
		ct->Column = ct;
		ct->Left = lastleft;
		lastleft->Right = ct;
		lastleft = ct;
		columns.push_back(ct);
		elements.push_back(ct);
	}
	lastleft->Right = &root;
	root.Left = lastleft;

	for (int i = 0; i < tvset.size(); i++) {
		vector<int> curtv = tvset[i];
		vector<DLX_column *> tvrow;
		for (int j = 0; j < curtv.size(); j++) {
			DLX_column *ctve;
			ctve = new (DLX_column);
			//column corresponds to characteristic vector of LS or smth of that kind
			int k = j*n + curtv[j];

			ctve->Column = columns[k];
			ctve->Column->size++;
			ctve->Down = columns[k];
			ctve->Up = columns[k]->Up;
			ctve->Up->Down = ctve;
			ctve->Down->Up = ctve;
			ctve->row_id = i;
			//	ctve->column_number = k;
			ctve->size = -10;
			elements.push_back(ctve);
			tvrow.push_back(ctve);
		}

		for (int j = 0; j < tvrow.size() - 1; j++) {
			tvrow[j]->Right = tvrow[j + 1];
			tvrow[j]->Right->Left = tvrow[j];
		}
		tvrow[tvrow.size() - 1]->Right = tvrow[0];
		tvrow[0]->Left = tvrow[tvrow.size() - 1];
	}
	DLX_column *pr = &root;

}

void DLX_orth::find_all_transversals(int k, DLX_column &h, vector<DLX_column*> &ps, vector<transversal_t> &tvr) {
	//pd = partial solution
	if (h.Right == &h) {
		vector<int> tmpv;
		for (int i = 0; i < ps.size(); i++) {
			tmpv.push_back(ps[i]->row_id);
		}
		tvr.push_back(tmpv);
		//cout << tvr.size() << endl;
		//print_solution(ps);
	}
	else {
		DLX_column * c = NULL;
		choose_c(h, c);
		//	cout << "picked column " << c->column_number << endl;
		cover(c);
		DLX_column * r = c->Down;
		while (r != c) {
			ps.push_back(r);
			DLX_column * j;
			j = r->Right;
			while (j != r) {
				cover(j->Column);
				j = j->Right;
			}

			find_all_transversals(k + 1, h, ps, tvr);

			r = ps.back();
			//questionable.
			ps.pop_back();
			c = r->Column;

			j = r->Left;
			while (j != r) {
				uncover(j->Column);
				j = j->Left;
			}
			r = r->Down;
		}
		uncover(c);
	}
}

vector<vector<int>> DLX_orth::find_tv_dlx(const latinsquare_t square, const bool is_diag) {
	assert(is_latinsquare(square));
	const int n = square.size();
	
	DLX_column *root;
	root = new (DLX_column);
	vector<DLX_column*> elements;
	square_to_DLX(*root, square, elements, is_diag);
	vector<DLX_column*> ps;
	ps.clear();
	vector<vector<int>> tvr;
	find_all_transversals(0, *root, ps, tvr);
	
	//cout << "Found " << tvr.size() << " transversals\n";
	for (int i = 0; i < tvr.size(); i++) {		
		sort(tvr[i].begin(), tvr[i].end());
		for (int j = 0; j < tvr[i].size(); j++) {
			tvr[i][j] = tvr[i][j] % n;
		}
		//printvector(tvr[i]);
		//cout << endl;
	}
	for (auto i = 0; i < elements.size(); i++) {
		delete elements[i];
	}
	elements.clear();
	return tvr;
}

// Find all orthogonal mates for a given Latin square:
vector<latinsquare_t> DLX_orth::find_all_orth_mates(const latinsquare_t square) {
	const unsigned n = square.size();
	assert(is_latinsquare(square));
	bool is_diag = true;
	vector<vector<int>> trm = find_tv_dlx(square, is_diag);
	DLX_column *root;
	root = new (DLX_column);
	vector<DLX_column*> elements;
	transversals_to_dlx(*root, trm, elements);
	vector<DLX_column*> ps;
	ps.clear();
	vector<transversal_t> transversals;
	find_all_transversals(0, *root, ps, transversals);
	for (int i = 0; i < transversals.size(); i++) {
		sort(transversals[i].begin(), transversals[i].end());
	}

	for (int i = 0; i < elements.size(); i++) {
		delete elements[i];
	}

	vector<latinsquare_t> orth_mates;
	if (transversals.size() > 0) {
		//out << transversals.size() << " sets of disjoint transversals" << endl;

		for (int i = 0; i < transversals.size(); i++) {
			latinsquare_t orth_square(n, row_t(n));
			for (unsigned u = 0; u < n; u++) {
				for (unsigned v = 0; v < n; v++) {
					orth_square[v][trm[transversals[i][u]][v]] = u;
				}
			}
			assert(is_latinsquare(orth_square));
			orth_mates.push_back(orth_square);
		}

	}

	return orth_mates;
}

// Find all orthogonal mates for a given Latin square:
LS_result DLX_orth::find_transversals_and_orth_mates(const latinsquare_t square) {
	const unsigned n = square.size();
	assert(is_latinsquare(square));
	// Find all diagonal transversals:
	vector<vector<int>> diag_transversals = find_tv_dlx(square, true);
	// Find all transversals:
	vector<vector<int>> transversals = find_tv_dlx(square, false);
	DLX_column *root;
	root = new (DLX_column);
	vector<DLX_column*> elements;
	transversals_to_dlx(*root, diag_transversals, elements);
	vector<DLX_column*> ps;
	ps.clear();
	// Find disjoint sets of transversals:
	vector<transversal_t> disjoint_transversals_sets;
	find_all_transversals(0, *root, ps, disjoint_transversals_sets);
	for (int i = 0; i < disjoint_transversals_sets.size(); i++) {
		sort(disjoint_transversals_sets[i].begin(), disjoint_transversals_sets[i].end());
	}

	for (int i = 0; i < elements.size(); i++) {
		delete elements[i];
	}

	vector<latinsquare_t> orth_mates;
	if (disjoint_transversals_sets.size() > 0) {
		//out << disjoint_transversals_sets.size() << " sets of disjoint transversals" << endl;

		for (int i = 0; i < disjoint_transversals_sets.size(); i++) {
			latinsquare_t orth_square(n, row_t(n));
			for (unsigned u = 0; u < n; u++) {
				for (unsigned v = 0; v < n; v++) {
					orth_square[v][diag_transversals[disjoint_transversals_sets[i][u]][v]] = u;
				}
			}
			assert(is_latinsquare(orth_square));
			orth_mates.push_back(orth_square);
		}

	}
	LS_result ls_res;
	ls_res.orth_mates = orth_mates;
	ls_res.diag_transv = diag_transversals.size();
	ls_res.transv = transversals.size();

	return ls_res;
}
