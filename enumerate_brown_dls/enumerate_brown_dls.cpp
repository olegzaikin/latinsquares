#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <set>
#include <sstream>

#include "../dlx_mols/dlx_orth.h"

using namespace std;

#define row_t vector<int>
#define matrix_t vector<row_t>

string program = "enumerate_brown_dls";
string version = "0.0.1";

bool is_latin_rows(const row_t row0, const row_t row1) {
    const unsigned n = row0.size();
    assert(n > 0 and n < 11);
    assert(row1.size() == n);
    bool res = false;
    for (unsigned i=0; i<n; i++) {
        if (row0[i] == row1[i]) return false; 
    }
    return true;
}

bool is_latin_square(const matrix_t ls) {
    const unsigned n = ls.size();
    assert(n > 0 and n < 11);
    bool res = true;
    const vector<unsigned> ones_vec(n, 1);
    for (unsigned i=0; i<n; i++) {
        vector<unsigned> values_num_vec(n, 0);
        for (unsigned j=0; j<n; j++) {
            values_num_vec[ls[i][j]]++;
        }
        if (values_num_vec != ones_vec) return false;
    }
    for (unsigned i=0; i<n; i++) {
        vector<unsigned> values_num_vec(n, 0);
        for (unsigned j=0; j<n; j++) {
            values_num_vec[ls[j][i]]++;
        }
        if (values_num_vec != ones_vec) return false;
    }
    return res;
}

bool is_diag_latin_square(const matrix_t ls) {
    const unsigned n = ls.size();
    assert(n > 0 and n < 11);
    bool res = true;
    if (not is_latin_square(ls)) return false;
    const vector<unsigned> ones_vec(n, 1);
    vector<unsigned> values_num_main_diag(n, 0);
    for (unsigned i=0; i<n; i++) {
        values_num_main_diag[ls[i][i]]++;
    }
    if (values_num_main_diag != ones_vec) return false;
    vector<unsigned> values_num_main_antidiag(n, 0);
    for (unsigned i=0; i<n; i++) {
        values_num_main_antidiag[ls[n-1-i][i]]++;
    }
    if (values_num_main_antidiag != ones_vec) return false;
    return res;
}

void print(matrix_t matrix) {
    for (auto &row : matrix) {
        for (auto &x : row) cout << x << " ";
        cout << endl;
    }
    cout << endl;
}

matrix_t turn_square_brown_style(const matrix_t quarter_ls) {
    const unsigned quarter_ls_size = quarter_ls.size();
    const unsigned n = quarter_ls_size * 2;
    assert(n > 0 and n < 11);
    matrix_t ls(n, row_t(n, -1));
    for (unsigned i=0; i<quarter_ls_size; i++) {
        // A top-left:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i][j] = quarter_ls[i][j];
        }
        // B = 9-A, bottom-left:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i+quarter_ls_size][j] = n - 1 - ls[i][j];
        }
        // B top-right:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i][j+quarter_ls_size] = n - 1 - ls[i][j];
        }
        // A bottom-right:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i+quarter_ls_size][j+quarter_ls_size] = ls[i][j];
        }
    }
    assert(is_latin_square(ls));
    return ls;
}

matrix_t turn_square_dmd_style(const matrix_t quarter_ls) {
    const unsigned quarter_ls_size = quarter_ls.size();
    const unsigned n = quarter_ls_size * 2;
    assert(n > 0 and n < 11);
    matrix_t ls(n, row_t(n, -1));
    for (unsigned i=0; i<quarter_ls_size; i++) {
        // A top-left:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i][j] = quarter_ls[i][j];
        }
        // B = 9-A, bottom-left:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i+quarter_ls_size][j] = n - 1 - ls[i][j];
        }
        // B' = reverse_columns(B), top-right:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i][j+quarter_ls_size] = n - 1 - ls[i][quarter_ls_size-1-j];
        }
        // A' = reverse_columns(A), bottom-right:
        for (unsigned j=0; j<quarter_ls_size; j++) {
            ls[i+quarter_ls_size][j+quarter_ls_size] = ls[i][quarter_ls_size-1-j];
        }
    }
    assert(is_latin_square(ls));
    return ls;
}

// Algorithm from 'On the Construction of Triples of Diagonal Latin Squares of Order 10':
vector<matrix_t> generate_ls_order5() {
    const unsigned n = 5;
    row_t perm;
    for (short int i=0; i<n; i++) perm.push_back(i);
    vector<row_t> perumtations;
    do {
       perumtations.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));
    cout << perumtations.size() << " permutations of size " << n << endl;
    vector<matrix_t> ls_arr;
    matrix_t ls(n, row_t(n, -1));
    for (auto &row0 : perumtations) {
        for (auto &row1 : perumtations) {
            if (not is_latin_rows(row0,row1)) continue;
            for (auto &row2 : perumtations) {
                if ((not is_latin_rows(row0,row2)) or 
                    (not is_latin_rows(row1,row2))) continue;
                for (auto &row3 : perumtations) {
                    if ((not is_latin_rows(row0,row3)) or 
                        (not is_latin_rows(row1,row3)) or 
                        (not is_latin_rows(row2,row3))) continue;
                    for (auto &row4 : perumtations) {
                        if ((not is_latin_rows(row0,row4)) or 
                            (not is_latin_rows(row1,row4)) or 
                            (not is_latin_rows(row2,row4)) or 
                            (not is_latin_rows(row3,row4))) continue;
                        ls[0] = row0;
                        ls[1] = row1;
                        ls[2] = row2;
                        ls[3] = row3;
                        ls[4] = row4;
                        if (is_latin_square(ls)) ls_arr.push_back(ls);
                    }
                }
            }
        }
    }
    return ls_arr;
}

// Normalize DLS, i.e. make its main diagonal 0, 1, ..., n-1
matrix_t normalize_main_diag(const matrix_t mtrx) {
    const unsigned n = mtrx.size();
    assert(n > 0 and n < 11);
    matrix_t norm_mtrx = mtrx;
    vector<unsigned> norm_perm(n, 0);
    for (unsigned i = 0; i < n; i++) norm_perm[mtrx[i][i]] = i;
    for (auto &row : norm_mtrx) {
        for (unsigned i = 0; i < n; i++) {
            if (row[i] == -1) continue;
            row[i] = norm_perm[row[i]];
        }
    }
    //cout << "*" << endl;
    //print(mtrx);
    //print(norm_mtrx);
    for (unsigned i = 0; i < n; i++) {
        assert(norm_mtrx[i][i] == i);
    }
    return norm_mtrx;
}

// Calculate the orhtogonality characteristics of a pair of Latin squares:
unsigned calc_orth_char(const latinsquare_t dls1, const latinsquare_t dls2) {
	unsigned orth_char = 0;
	std::set<string> ordered_pairs;
	size_t n = dls1.size();
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			stringstream sstream;
			sstream << dls1[i][j] << dls2[i][j];
			ordered_pairs.insert(sstream.str());
		}
	}
	return ordered_pairs.size();
}

int main(int argc, char *argv[])
{
    vector<string> argv_str;
    for (unsigned i=0; i<argc; i++) argv_str.push_back(argv[i]);
	if ((argc == 2) and (argv_str[1] == "-v")) {
        cout << program << " of version " << version << endl;
        return 1;
    }

    cout << "Running " << program << " of version " << version << endl;

    vector<matrix_t> ls_arr = generate_ls_order5();
    cout << ls_arr.size() << " LS of order 5" << endl;

    const unsigned n = 10;
    row_t perm;
    for (short int i=0; i<n; i++) perm.push_back(i);
    vector<row_t> perumtations;
    do {
       perumtations.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));
    cout << perumtations.size() << " permutations of size " << n << endl;

    matrix_t quarter_ls;
    // A 5x5 square from DMD 2016: 
    quarter_ls = {{0, 1, 2, 3, 4},
                  {1, 2, 3, 4, 9},
                  {2, 3, 4, 0, 8},
                  {3, 5, 9, 8, 2},
                  {4, 0, 8, 7, 6}};

    //matrix_t ls = turn_square_brown_style(quarter_ls);
    matrix_t ls = turn_square_dmd_style(quarter_ls);
    cout << "Turn LS :" << endl;
    print(ls);
    unsigned dls_num = 0;
    unsigned k=0;
    cout << "Making DLS by permuting rows and columns" << endl;
    vector<matrix_t> dls_arr;
    for (auto &perm : perumtations) {
        matrix_t perm_row_ls(n, row_t(n,-1));
        matrix_t perm_col_ls(n, row_t(n,-1));
        for (unsigned i=0; i<n; i++) {
            perm_row_ls[i] = ls[perm[i]];
        }
        assert(is_latin_square(perm_row_ls));
        if (is_diag_latin_square(perm_row_ls)) {
            dls_num++;
            dls_arr.push_back(perm_row_ls);
        }
        for (unsigned i=0; i<n; i++) {
            for (unsigned j=0; j<n; j++) {
                perm_col_ls[i][j] = ls[i][perm[j]];
            }
        }
        assert(is_latin_square(perm_col_ls));
        if (is_diag_latin_square(perm_col_ls)) {
            dls_num++;
            dls_arr.push_back(perm_col_ls);
        }
        k++;
        if (k % 1000000 == 0) cout << k << " processed" << endl;
    }
    cout << k << " processed" << endl;

    cout << dls_arr.size() << " DLS" << endl;

    for (auto &dls : dls_arr) assert(DLX_orth::is_diag_latinsquare(dls));

    set<matrix_t> normalized_dls_set;
    for (auto &dls : dls_arr) {
        //print(dls);
        matrix_t norm_dls = normalize_main_diag(dls);
        normalized_dls_set.insert(norm_dls);
    }
    cout << normalized_dls_set.size() << " normalized DLS" << endl;
    dls_arr.clear();

	unsigned max_orth_char = 0;
    unsigned max_orth_triples_num = 0;
    k=0;
    unsigned no_orth_mate_dls_num = 0;
    unsigned one_orth_mate_dls_num = 0;
    unsigned two_orth_mate_dls_num = 0;
    unsigned three_orth_mate_dls_num = 0;
	for (auto &norm_dls : normalized_dls_set) {
		vector<latinsquare_t> orth_mates = DLX_orth::find_all_orth_mates(norm_dls);
        k++;
        if (k % 1000 == 0) cout << k << " normalized DLS processed" << endl;
		// For all pairs of DLS which are orthogonal to the current square and form a triple:
        if (orth_mates.size() == 0) no_orth_mate_dls_num++;
        if (orth_mates.size() == 1) one_orth_mate_dls_num++;
        if (orth_mates.size() == 2) two_orth_mate_dls_num++;
        if (orth_mates.size() == 3) three_orth_mate_dls_num++;
		if (orth_mates.size() < 2) continue;
		for (unsigned j = 0; j < orth_mates.size() - 1; j++) {
			for (unsigned j2 = j+1; j2 < orth_mates.size(); j2++) {
				unsigned orth_char = calc_orth_char(orth_mates[j], orth_mates[j2]);
				//if (orth_char == max_orth_char) cout << orth_char << endl;
                if (orth_char == 74) {
                    max_orth_triples_num++;
                    cout << "max_orth_triples_num : " << max_orth_triples_num << endl;
                }
				if (max_orth_char == 0 or orth_char > max_orth_char) {
					max_orth_char = orth_char;
					cout << "Updated max_orth_char : " << max_orth_char << endl;
					print(norm_dls);
				}
			}
		}
	}
    cout << k << " normalized DLS processed" << endl;
    cout << no_orth_mate_dls_num    << " DLS with 0 orthogonal mates" << endl;
    cout << one_orth_mate_dls_num   << " DLS with 1 orthogonal mates" << endl;
    cout << two_orth_mate_dls_num   << " DLS with 2 orthogonal mates" << endl;
    cout << three_orth_mate_dls_num << " DLS with 3 orthogonal mates" << endl;

    return 0;
}
