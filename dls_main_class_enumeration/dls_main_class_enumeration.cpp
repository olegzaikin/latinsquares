// Created on: 23 Oct 2024
// Author: Oleg Zaikin
// E-mail: oleg.zaikin@icc.ru
//
// For a given order n and a file of ESODLS CMS, generate all ESODLS CMS of
// order n and compare them with that from the file.
// 
// Example:
//   ./dls_main_class_enumeration 7 ./cms_n7
//==========================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <bitset>
#include <set>

using namespace std;

#define row_t vector<short int>
#define matrix_t vector<row_t>

string prog = "dls_main_class_enumeration";
string version = "0.1.1";

void print(matrix_t matrix) {
    for (auto &row : matrix) {
        for (auto &x : row) {
            if (x >= 0 and x < 10) cout << "0";
            cout << x << " ";
        }
        cout << endl;
    }
    cout << endl;
}

bool is_digits(const string str)
{
    return all_of(str.begin(), str.end(), ::isdigit);
}

short int strtoi(string s) {
	assert(not s.empty());
	short int x = atoi(s.c_str());
	return x;
}

bool is_diag_cms(const matrix_t cms) {
    const unsigned n = cms.size();
    assert(n > 0 and n < 11); 
    vector<unsigned> cms_values_num(n*n, 0);
    // Make the correct array of two diagonals:
    row_t sorted_two_diag_dls_cms(n*2, -1);
    for (unsigned i=0; i<n; i++) {
        sorted_two_diag_dls_cms[i] = i*n + i;
        sorted_two_diag_dls_cms[n + i] = (n-1-i)*n + i;
    }
    sort(sorted_two_diag_dls_cms.begin(), sorted_two_diag_dls_cms.end());
    row_t two_diag_dls_cms(2*n, -1);
    for (unsigned i=0; i<n; i++) {
        two_diag_dls_cms[i] = cms[i][i];
        two_diag_dls_cms[n + i] = cms[n-1-i][i];
        for (unsigned j=0; j<n; j++) {
            assert(cms[i][j] < cms_values_num.size());
            cms_values_num[cms[i][j]]++;
        }
    }
    sort(two_diag_dls_cms.begin(), two_diag_dls_cms.end());
    if (two_diag_dls_cms != sorted_two_diag_dls_cms) return false;
    for (unsigned i=0; i<n; i++) {
        if (cms_values_num[i] == 0) return false;
    }
    return true;
}

// Read CMSs from a given file:
set<matrix_t> read_cms(const string filename, const unsigned n){
	set<matrix_t> cms_set;
	ifstream in;
	in.open(filename);
	assert(in.is_open());
	string s;
	matrix_t cms;
	while (getline(in, s)) {
		if (s == "" and cms.size() > 0){
			assert(is_diag_cms(cms));
			cms_set.insert(cms);
			cms.clear();
		}
		else {
			row_t tmp;
            stringstream sstream(s);
            string word;
            while (sstream >> word) {
                if (not is_digits(word)) {
                    tmp.clear();
                    break;
                }
				tmp.push_back(strtoi(word));
            }
            if (not tmp.empty()) cms.push_back(tmp);
		}
	}
	in.close();
	return cms_set;
}

vector<matrix_t> read_dls_string(const string fname, const unsigned n) {
    assert(fname != "");
    vector<matrix_t> diag_ls_lst;

    ifstream f(fname);
    string str;
    while (getline(f, str)) {
        str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
        str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
        if (str.size() != n*n) continue;
        //cout << str << endl;
        matrix_t dls;
        row_t row;
        for (unsigned i=0; i<str.size(); i++) {
            row.push_back(std::stoi(str.substr(i,1))); // row.push_back(std::stoi(str[i]));
            if (row.size() == n) {
                dls.push_back(row);
                row.clear();
            }
        }
        diag_ls_lst.push_back(dls);
    }
    f.close();
    return diag_ls_lst;
}

// Normalize matrix (LS or CMS), i.e. make its main diagonal 0, 1, ..., n-1
matrix_t normalize(const matrix_t mtrx) {
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

void swap_rows(matrix_t &dls, const unsigned i1, const unsigned i2) {
    const unsigned n = dls.size();
    assert(n > 0 and n < 11);
    assert(i1 < n and i2 < n);
    swap(dls[i1], dls[i2]);
}

void swap_columns(matrix_t &dls, const unsigned j1, const unsigned j2) {
     const unsigned n = dls.size();
    assert(n > 0 and n < 11);
    assert(j1 < n and j1 < n);
    for (unsigned i = 0; i < n; i++) swap(dls[i][j1], dls[i][j2]);
}

set<matrix_t> rotate(const matrix_t dls) {
    const unsigned n = dls.size();

    assert(n > 0 and n < 11); 
    set<matrix_t> rotate_dls_set;

    // 90 degree rotation:
    matrix_t one_rotate_dls = dls;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            one_rotate_dls[j][n-1-i] = dls[i][j];
        }
    }

    matrix_t two_rotate_dls = one_rotate_dls;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            two_rotate_dls[j][n-1-i] = one_rotate_dls[i][j];
        }
    }

    matrix_t three_rotate_dls = two_rotate_dls;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            three_rotate_dls[j][n-1-i] = two_rotate_dls[i][j];
        }
    }

    assert(is_diag_cms(one_rotate_dls));
    assert(is_diag_cms(two_rotate_dls));
    assert(is_diag_cms(three_rotate_dls));

    rotate_dls_set.insert(dls);
    rotate_dls_set.insert(one_rotate_dls);
    rotate_dls_set.insert(two_rotate_dls);
    rotate_dls_set.insert(three_rotate_dls);

    assert(rotate_dls_set.size() == 4);

    return rotate_dls_set;
}

set<matrix_t> reflect(const matrix_t dls) {
    const unsigned n = dls.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> reflect_dls_set;

    // Reflect horizontally across the y-axis:
    matrix_t reflect_horiz_dls = dls;
    unsigned k1 = (unsigned)floor((double)n/2);
    unsigned k2 = (unsigned)ceil((double)n/2);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < k1; j++) {
            reflect_horiz_dls[i][j] = dls[i][n-j-1];
        }
        for (unsigned j = k2; j < n; j++) {
            reflect_horiz_dls[i][j] = dls[i][k1 - 1 - (j-k2)];
        }
    }

    // Reflect vertically across the x-axis:
    matrix_t reflect_vert_dls = dls;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < k1; j++) {
            reflect_vert_dls[j][i] = dls[n-j-1][i];
        }
        for (unsigned j = k2; j < n; j++) {
            reflect_vert_dls[j][i] = dls[k1 - 1 - (j-k2)][i];
        }
    }

    // Reflect across the main diagonal == transposition:
    matrix_t reflect_maindiag_dls = dls;
    for (unsigned i = 1; i < n; i++) {
        for (unsigned j = 0; j < i; j++) {
            swap(reflect_maindiag_dls[i][j], reflect_maindiag_dls[j][i]);
        }
    }

    // Reflect across the main antidiagonal:
    matrix_t reflect_antidiag_dls = dls;
    for (unsigned i = 0; i < n-1; i++) {
        for (unsigned j = 0; j < n-1-i; j++) {
            swap(reflect_antidiag_dls[i][j], reflect_antidiag_dls[n-1-j][n-1-i]);
        }
    }

    //cout << endl;
    //print(dls);
    //print(reflect_horiz_dls);
    assert(is_diag_cms(reflect_horiz_dls));
    assert(is_diag_cms(reflect_vert_dls));
    assert(is_diag_cms(reflect_maindiag_dls));
    assert(is_diag_cms(reflect_antidiag_dls));

    reflect_dls_set.insert(dls);
    reflect_dls_set.insert(reflect_horiz_dls);
    reflect_dls_set.insert(reflect_vert_dls);
    reflect_dls_set.insert(reflect_maindiag_dls);
    reflect_dls_set.insert(reflect_antidiag_dls);

    assert(reflect_dls_set.size() == 5);

    return reflect_dls_set;
}

set<matrix_t> reflect_rotate(const matrix_t cms) {
    set<matrix_t> refl_rot_cms_set;
    set<matrix_t> reflection_cms_set = reflect(cms);
    for (auto &refl_cms : reflection_cms_set) {
        set<matrix_t> rotation_cms_set = rotate(refl_cms);
        for (auto &rot_cms : rotation_cms_set) {
            assert(is_diag_cms(rot_cms));
            refl_rot_cms_set.insert(rot_cms);
        }
    }
    return refl_rot_cms_set;
}

// Swap two rows symmetrically and two columns with the same indices.
// There are 2^(floor(n/2)) such variants where n is the order of a given
// diagonal Latin square.
set<matrix_t> tworows_twocols_symm(const matrix_t dls) {
    const unsigned n = dls.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> tworows_twocols_symm_dls_set;
    const unsigned subset_num = (unsigned)pow(2,n/2);
    for (unsigned i=0; i < subset_num; i++) {
        bitset<5> b{i};
        string str_bit_repres = b.to_string();
        str_bit_repres = str_bit_repres.substr(5 - n/2, n/2);
        matrix_t symm_dls = dls;
        for (unsigned row_ind = 0; row_ind < str_bit_repres.size(); row_ind++) {
            assert(row_ind < n/2);
            if (str_bit_repres[row_ind] == '1') {
                // Swap two rows symmetrically:
                swap_rows(symm_dls, row_ind, n-1-row_ind);
                // Swap two columns with the same indices:
                swap_columns(symm_dls, row_ind, n-1-row_ind);
            }
        }
        assert(is_diag_cms(symm_dls));
        tworows_twocols_symm_dls_set.insert(symm_dls);
    }
    assert(tworows_twocols_symm_dls_set.size() == subset_num);
    return tworows_twocols_symm_dls_set;
}

// Swap four rows symmetrically (2 + 2) and four columns with the same indices.
// There are (floor(n/2))! such variants where n is the order of a given
// diagonal Latin square.
set<matrix_t> fourrows_fourcols_symm(const matrix_t dls) {
    const unsigned n = dls.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> fourrows_fourcols_symm_dls_set;
    vector<unsigned> cur_indices;
    for (unsigned i=0; i<n/2; i++) cur_indices.push_back(i);
    unsigned perm_num = 0;
    do {
        matrix_t symm_dls_temp = dls;
        // Assign rows:
        for (unsigned i=0; i < cur_indices.size(); i++) {
            symm_dls_temp[i] = dls[cur_indices[i]];
            symm_dls_temp[n-1-i] = dls[n-1-cur_indices[i]];
        }
        // Assign columns:
        matrix_t symm_dls = symm_dls_temp;
        for (unsigned i=0; i < cur_indices.size(); i++) {
            for (unsigned j = 0; j < n; j++) {
                symm_dls[j][i] = symm_dls_temp[j][cur_indices[i]];
                symm_dls[j][n-1-i] = symm_dls_temp[j][n-1-cur_indices[i]];
            }
        }
        assert(is_diag_cms(symm_dls));
        fourrows_fourcols_symm_dls_set.insert(symm_dls);
        perm_num++;
    } while (std::next_permutation(cur_indices.begin(), cur_indices.end()));
    assert(fourrows_fourcols_symm_dls_set.size() == perm_num);
    return fourrows_fourcols_symm_dls_set;
}

set<matrix_t> reflect_rotate_symm(const matrix_t cms) {
    const unsigned n = cms.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> cms_set;
    set<matrix_t> tworows_twocols_symm_cms_set = tworows_twocols_symm(cms);
    for (auto &two_symm_cms : tworows_twocols_symm_cms_set) {
        set<matrix_t> fourrows_fourcols_symm_cms_set = fourrows_fourcols_symm(two_symm_cms);
        for (auto &four_symm_cms : fourrows_fourcols_symm_cms_set) {
            set<matrix_t> refl_rot_cms_set = reflect_rotate(four_symm_cms);
            for (auto &refl_rot_cms : refl_rot_cms_set) {
                assert(is_diag_cms(refl_rot_cms));
                cms_set.insert(refl_rot_cms);
            }
        }
    }
    return cms_set;
}

vector<matrix_t> generate_x_based_basic_partial_dls(const unsigned n) {
    assert(n > 0 and n < 11);
    vector<matrix_t> x_based_basic_partial_dls;
    vector<unsigned> main_diag;
    for (unsigned i=0; i<n; i++) main_diag.push_back(i);
    vector<unsigned> main_antidiag;
    for (unsigned i=0; i<n; i++) main_antidiag.push_back(i);
    do {
        bool is_diag = true;
        if (n % 2 == 1) {
            if (main_diag[n/2] != main_antidiag[n/2]) {
                is_diag = false;
                continue;
            }
            for (unsigned i=0; i<n; i++) {
                if (i == n/2) continue;
                if (main_diag[i] == main_antidiag[i] or main_diag[i] == main_antidiag[n-1-i]) {
                    is_diag = false;
                    continue;
                }
            }
        }
        else {
            for (unsigned i=0; i<n; i++) {
                if (main_diag[i] == main_antidiag[i] or main_diag[i] == main_antidiag[n-1-i]) {
                    is_diag = false;
                    continue;
                }
            }
        }
        if (is_diag) {
            matrix_t x_based_partial_dls(n,row_t(n,-1));
            for (unsigned i=0; i<n; i++) {
                x_based_partial_dls[i][i] = main_diag[i];
                x_based_partial_dls[n-1-i][i] = main_antidiag[i];
            }
            x_based_basic_partial_dls.push_back(x_based_partial_dls);
        }
    } while (std::next_permutation(main_antidiag.begin(), main_antidiag.end()));
    return x_based_basic_partial_dls;
}

matrix_t apply_all_cms(matrix_t x_based_partial_dls, set<matrix_t> esodls_cms_set) {
    const unsigned n = x_based_partial_dls.size();
    assert(n > 0 and n < 11);
    matrix_t new_x_based_partial_dls(n, row_t(n,-1));
    matrix_t norm_new_x_based_partial_dls(n, row_t(n,-1));
    matrix_t minimal_main_class_repres;
    vector<short int> minimal_two_diags_rows_by_rows;
    //cout << "basic x :" << endl;
    //print(x_based_partial_dls);
    for (auto cms : esodls_cms_set) {
        for (unsigned i=0; i<n; i++) {
            unsigned i2 = cms[i][i] / n;
            unsigned j2 = cms[i][i] % n;
            assert(i2 == j2 or i2 == n-1-j2 or j2 == n-1-i2);
            new_x_based_partial_dls[i2][j2] = x_based_partial_dls[i][i];
            i2 = cms[n-1-i][i] / n;
            j2 = cms[n-1-i][i] % n;
            assert(i2 == j2 or i2 == n-1-j2 or j2 == n-1-i2);
            new_x_based_partial_dls[i2][j2] = x_based_partial_dls[n-1-i][i];
        }
        // Normalize:
        norm_new_x_based_partial_dls = normalize(new_x_based_partial_dls);
        // Get the antidiag:
        vector<short int> two_diags_rows_by_rows(n*2, -1);
        unsigned k=0;
        for (unsigned i=0; i<n; i++) {
            for (unsigned j=0; j<n; j++) {
                if (norm_new_x_based_partial_dls[i][j] == -1) continue;
                assert(k < two_diags_rows_by_rows.size());
                two_diags_rows_by_rows[k++] = norm_new_x_based_partial_dls[i][j];
            }
        }
        two_diags_rows_by_rows.resize(k);
        //cout << "*" << endl;
        //print(norm_new_x_based_partial_dls);
        //for (auto &x : two_diags_rows_by_rows) cout << x << " ";
        //cout << endl;
        
        if (minimal_two_diags_rows_by_rows.empty() or two_diags_rows_by_rows < minimal_two_diags_rows_by_rows) {
            minimal_main_class_repres = norm_new_x_based_partial_dls;
            minimal_two_diags_rows_by_rows = two_diags_rows_by_rows;
        }
    }
    return minimal_main_class_repres;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cout << "Usage : " << prog << " DLS-order CMS-file" << endl;
        return 1;
    }
    const unsigned n = atoi(argv[1]);
    string cms_fname = argv[2];

    cout << prog << " of version " << version << " is running" << endl;
    cout << "DLS order : " << n << endl;
    cout << "CMS file name : " << cms_fname << endl;
    assert(n > 0 and n < 11);

    set<matrix_t> cms_from_file = read_cms(cms_fname, n);
    cout << cms_from_file.size() << " cms were read" << endl;
    //for (auto &cms : cms_from_file) print(cms);

    matrix_t trivial_cms(n,row_t(n, 0));
    for (unsigned i=0; i<n; i++) {
        for (unsigned j=0; j<n; j++) {
            trivial_cms[i][j] = i*n + j;
        }
    }

    //vector<matrix_t> dls_arr = read_dls_string(fname, n);

    cout << "Trivial CMS :" << endl;
    print(trivial_cms);

    //set<matrix_t> refl_rot_cms_set = reflect_rotate(trivial_cms);
    //cout << "refl_rot_cms_set size : " << refl_rot_cms_set.size() << endl;

    set<matrix_t> esodls_cms_set = reflect_rotate_symm(trivial_cms); 
    cout << "esodls_cms_set size : " << esodls_cms_set.size() << endl;

    set<matrix_t> diff_set;
    std::set_difference(
        cms_from_file.begin(), cms_from_file.end(),
        esodls_cms_set.begin(), esodls_cms_set.end(),
        std::inserter(diff_set, diff_set.end()));

    assert(cms_from_file.size() == esodls_cms_set.size() and diff_set.empty());

    cout << "CMS set is correct" << endl;

    stringstream sstream_cms;
    sstream_cms << "cms_esodls_n" << n << "_new.txt";
    string out_cms_filename = sstream_cms.str();
    cout << "Wriring ESODLS CMS to file " << out_cms_filename << endl;
    ofstream out_cms_file(out_cms_filename, ios_base::out);
    out_cms_file << "Writing " << esodls_cms_set.size() << " ESODLS CMS of order " << n << endl;
    unsigned k = 0;
    for (auto &cms : esodls_cms_set) {
        out_cms_file << "# " << k << endl;
        for (auto &row : cms) {
            for (unsigned j=0; j<row.size(); j++) {
                assert(row[j] >= 0 and row[j] < n*n);
                if (row[j] < 10) out_cms_file << " ";
                out_cms_file << row[j];
                if (j < row.size() - 1) out_cms_file << " ";
            }
            out_cms_file << endl;
        }
        out_cms_file << endl;
        k++;
    }
    out_cms_file.close();

    vector<matrix_t> x_based_basic_partial_dls = generate_x_based_basic_partial_dls(n);
    cout << x_based_basic_partial_dls.size() << " basic X-based fillings" << endl;

    set<matrix_t> main_class_repres_set;
    k=0;
    unsigned max_main_class_repres_set_size = 0;
    for (auto &x_partial_dls : x_based_basic_partial_dls) {
        matrix_t main_class_repres = apply_all_cms(x_partial_dls, esodls_cms_set);
        main_class_repres_set.insert(main_class_repres);
        if (main_class_repres_set.size() > max_main_class_repres_set_size) {
            max_main_class_repres_set_size = main_class_repres_set.size();
            cout << "main_class_repres_set size : " << main_class_repres_set.size() << endl;
        }
        k++;
        if (k % 10000 == 0) cout << k << " processed" << endl;
    }

    cout << "main_class_repres_set size : " << main_class_repres_set.size() << endl;

    stringstream sstream;
    sstream << "x_fillings_n" << n << "_new.txt";
    string x_fname = sstream.str();
    cout << "Wriring X-based partial fillings to file " << x_fname << endl;
    ofstream ofile(x_fname, ios_base::out);
    k = 0;
    for (auto &x_fill : main_class_repres_set) {
        ofile << "X-based filling " << k << ":\n";
        for (unsigned i=0; i < n; i++) {
            for (unsigned j=0; j < n; j++) {
                if (x_fill[i][j] == -1) ofile << "-";
                else ofile << x_fill[i][j];
                if (j < n-1) ofile << " ";
            }
            ofile << "\n";
        }
        ofile << "\n";
        k++;
    }
    ofile.close();

    return 0;
}
