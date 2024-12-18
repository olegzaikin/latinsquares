// Created on: 23 Oct 2024
// Author: Oleg Zaikin
// E-mail: oleg.zaikin@icc.ru
//
// For a given order n, a file with all DLS normalized by the main diagonal,
// and a file of all ESODLS CMS, generate all main classes for DLS, and then
// generate all ESODLS CMS of order n and compare them with that from the file.
// 
// Example:
//   ./dls_main_class_enumeration 7 ./dls_n7 ./cms_n7
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
#include <chrono>

using namespace std;

#define row_t vector<short int>
#define matrix_t vector<row_t>

string prog = "dls_main_class_enumeration";
string version = "0.2.9";

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
    const short int n = cms.size();
    assert(n > 0 and n < 11); 
    vector<short int> cms_values_num(n*n, 0);
    // Make the correct array of two diagonals:
    row_t sorted_two_diag_dls_cms(n*2, -1);
    for (short int i=0; i<n; i++) {
        sorted_two_diag_dls_cms[i] = i*n + i;
        sorted_two_diag_dls_cms[n + i] = (n-1-i)*n + i;
    }
    sort(sorted_two_diag_dls_cms.begin(), sorted_two_diag_dls_cms.end());
    row_t two_diag_dls_cms(2*n, -1);
    for (short int i=0; i<n; i++) {
        two_diag_dls_cms[i] = cms[i][i];
        two_diag_dls_cms[n + i] = cms[n-1-i][i];
        for (short int j=0; j<n; j++) {
            assert(cms[i][j] < cms_values_num.size());
            cms_values_num[cms[i][j]]++;
        }
    }
    sort(two_diag_dls_cms.begin(), two_diag_dls_cms.end());
    if (two_diag_dls_cms != sorted_two_diag_dls_cms) return false;
    for (short int i=0; i<n; i++) {
        if (cms_values_num[i] == 0) return false;
    }
    return true;
}

// Read CMSs from a given file:
set<matrix_t> read_cms(const string filename, const short int n){
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

vector<matrix_t> read_dls_string(const string fname, const short int n) {
    assert(fname != "");
    vector<matrix_t> dls_set;

    ifstream f(fname);
    string str;
    while (getline(f, str)) {
        str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
        str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
        if (str.size() != n*n) continue;
        matrix_t dls;
        row_t row;
        for (short int i=0; i<str.size(); i++) {
            row.push_back(std::stoi(str.substr(i,1))); // row.push_back(std::stoi(str[i]));
            if (row.size() == n) {
                dls.push_back(row);
                row.clear();
            }
        }
        dls_set.push_back(dls);
    }
    f.close();
    return dls_set;
}

// Normalize matrix (LS or CMS), i.e. make its main diagonal 0, 1, ..., n-1
matrix_t normalize_main_diag(const matrix_t mtrx) {
    const short int n = mtrx.size();
    assert(n > 0 and n < 11);
    matrix_t norm_mtrx = mtrx;
    vector<short int> norm_perm(n, 0);
    for (short int i = 0; i < n; i++) norm_perm[mtrx[i][i]] = i;
    for (auto &row : norm_mtrx) {
        for (short int i = 0; i < n; i++) {
            if (row[i] == -1) continue;
            row[i] = norm_perm[row[i]];
        }
    }
    for (short int i = 0; i < n; i++) {
        assert(norm_mtrx[i][i] == i);
    }
    return norm_mtrx;
}

// Normalize matrix (LS or CMS), i.e. make its first row 0, 1, ..., n-1
matrix_t normalize_first_row(const matrix_t mtrx) {
    const short int n = mtrx.size();
    assert(n > 0 and n < 11);
    matrix_t norm_mtrx = mtrx;
    vector<short int> norm_perm(n, 0);
    for (short int i = 0; i < n; i++) norm_perm[mtrx[0][i]] = i;
    for (auto &row : norm_mtrx) {
        for (short int i = 0; i < n; i++) {
            if (row[i] == -1) continue;
            row[i] = norm_perm[row[i]];
        }
    }
    for (short int i = 0; i < n; i++) {
        assert(norm_mtrx[0][i] == i);
    }
    return norm_mtrx;
}

void swap_rows(matrix_t &dls, const short int i1, const short int i2) {
    const short int n = dls.size();
    assert(n > 0 and n < 11);
    assert(i1 < n and i2 < n);
    swap(dls[i1], dls[i2]);
}

void swap_columns(matrix_t &dls, const short int j1, const short int j2) {
     const short int n = dls.size();
    assert(n > 0 and n < 11);
    assert(j1 < n and j1 < n);
    for (short int i = 0; i < n; i++) swap(dls[i][j1], dls[i][j2]);
}

set<matrix_t> rotate(const matrix_t dls) {
    const short int n = dls.size();

    assert(n > 0 and n < 11); 
    set<matrix_t> rotate_dls_set;

    // 90 degree rotation:
    matrix_t one_rotate_dls = dls;
    for (short int i = 0; i < n; i++) {
        for (short int j = 0; j < n; j++) {
            one_rotate_dls[j][n-1-i] = dls[i][j];
        }
    }

    matrix_t two_rotate_dls = one_rotate_dls;
    for (short int i = 0; i < n; i++) {
        for (short int j = 0; j < n; j++) {
            two_rotate_dls[j][n-1-i] = one_rotate_dls[i][j];
        }
    }

    matrix_t three_rotate_dls = two_rotate_dls;
    for (short int i = 0; i < n; i++) {
        for (short int j = 0; j < n; j++) {
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
    const short int n = dls.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> reflect_dls_set;

    // Reflect horizontally across the y-axis:
    matrix_t reflect_horiz_dls = dls;
    short int k1 = (short int)floor((double)n/2);
    short int k2 = (short int)ceil((double)n/2);
    for (short int i = 0; i < n; i++) {
        for (short int j = 0; j < k1; j++) {
            reflect_horiz_dls[i][j] = dls[i][n-j-1];
        }
        for (short int j = k2; j < n; j++) {
            reflect_horiz_dls[i][j] = dls[i][k1 - 1 - (j-k2)];
        }
    }

    // Reflect vertically across the x-axis:
    matrix_t reflect_vert_dls = dls;
    for (short int i = 0; i < n; i++) {
        for (short int j = 0; j < k1; j++) {
            reflect_vert_dls[j][i] = dls[n-j-1][i];
        }
        for (short int j = k2; j < n; j++) {
            reflect_vert_dls[j][i] = dls[k1 - 1 - (j-k2)][i];
        }
    }

    // Reflect across the main diagonal == transposition:
    matrix_t reflect_maindiag_dls = dls;
    for (short int i = 1; i < n; i++) {
        for (short int j = 0; j < i; j++) {
            swap(reflect_maindiag_dls[i][j], reflect_maindiag_dls[j][i]);
        }
    }

    // Reflect across the main antidiagonal:
    matrix_t reflect_antidiag_dls = dls;
    for (short int i = 0; i < n-1; i++) {
        for (short int j = 0; j < n-1-i; j++) {
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
    const short int n = dls.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> tworows_twocols_symm_dls_set;
    const short int subset_num = (short int)pow(2,n/2);
    for (short int i=0; i < subset_num; i++) {
        bitset<5> b{(unsigned long long)i};
        string str_bit_repres = b.to_string();
        str_bit_repres = str_bit_repres.substr(5 - n/2, n/2);
        matrix_t symm_dls = dls;
        for (short int row_ind = 0; row_ind < str_bit_repres.size(); row_ind++) {
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
    const short int n = dls.size();
    assert(n > 0 and n < 11); 
    set<matrix_t> fourrows_fourcols_symm_dls_set;
    vector<short int> cur_indices;
    for (short int i=0; i<n/2; i++) cur_indices.push_back(i);
    unsigned perm_num = 0;
    do {
        matrix_t symm_dls_temp = dls;
        // Assign rows:
        for (short int i=0; i < cur_indices.size(); i++) {
            symm_dls_temp[i] = dls[cur_indices[i]];
            symm_dls_temp[n-1-i] = dls[n-1-cur_indices[i]];
        }
        // Assign columns:
        matrix_t symm_dls = symm_dls_temp;
        for (short int i=0; i < cur_indices.size(); i++) {
            for (short int j = 0; j < n; j++) {
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
    const short int n = cms.size();
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

vector<matrix_t> generate_x_based_basic_partial_dls(const short int n) {
    assert(n > 0 and n < 11);
    vector<matrix_t> x_based_basic_partial_dls;
    vector<short int> main_diag;
    for (short int i=0; i<n; i++) main_diag.push_back(i);
    vector<short int> main_antidiag;
    for (short int i=0; i<n; i++) main_antidiag.push_back(i);
    do {
        bool is_diag = true;
        if (n % 2 == 1) {
            if (main_diag[n/2] != main_antidiag[n/2]) {
                is_diag = false;
                continue;
            }
            for (short int i=0; i<n; i++) {
                if (i == n/2) continue;
                if (main_diag[i] == main_antidiag[i] or main_diag[i] == main_antidiag[n-1-i]) {
                    is_diag = false;
                    continue;
                }
            }
        }
        else {
            for (short int i=0; i<n; i++) {
                if (main_diag[i] == main_antidiag[i] or main_diag[i] == main_antidiag[n-1-i]) {
                    is_diag = false;
                    continue;
                }
            }
        }
        if (is_diag) {
            matrix_t x_based_partial_dls(n,row_t(n,-1));
            for (short int i=0; i<n; i++) {
                x_based_partial_dls[i][i] = main_diag[i];
                x_based_partial_dls[n-1-i][i] = main_antidiag[i];
            }
            x_based_basic_partial_dls.push_back(x_based_partial_dls);
        }
    } while (std::next_permutation(main_antidiag.begin(), main_antidiag.end()));
    return x_based_basic_partial_dls;
}

matrix_t find_main_class_repres_all_cms(const matrix_t partial_dls, set<matrix_t> cms_set) {
    const short int n = partial_dls.size();
    assert(n > 0 and n < 11);
    matrix_t new_partial_dls(n, row_t(n,-1));
    matrix_t norm_partial_dls(n, row_t(n,-1));
    matrix_t minimal_main_class_repres;
    vector<short int> minimal_one_dim_partial_dls_repres;
    //cout << "basic x :" << endl;
    //print(x_based_partial_dls);
    for (auto cms : cms_set) {
        for (short int i=0; i<n; i++) {
            for (short int j=0; j<n; j++) {
                const short int i2 = cms[i][j] / n;
                const short int j2 = cms[i][j] % n;
                new_partial_dls[i2][j2] = partial_dls[i][j];
            }
        }
        // Normalize:
        norm_partial_dls = normalize_main_diag(new_partial_dls);
        // Get the antidiag:
        vector<short int> one_dim_partial_dls_repres(n*n, -1);
        unsigned k=0;
        for (short int i=0; i<n; i++) {
            for (short int j=0; j<n; j++) {
                if (norm_partial_dls[i][j] == -1) continue;
                assert(k < one_dim_partial_dls_repres.size());
                one_dim_partial_dls_repres[k++] = norm_partial_dls[i][j];
            }
        }
        one_dim_partial_dls_repres.resize(k);

        if (minimal_one_dim_partial_dls_repres.empty() or one_dim_partial_dls_repres < minimal_one_dim_partial_dls_repres) {
            minimal_main_class_repres = norm_partial_dls;
            minimal_one_dim_partial_dls_repres = one_dim_partial_dls_repres;
        }
    }
    return minimal_main_class_repres;
}

vector<matrix_t> find_all_dls_with_ascending_row_from_main_class(const matrix_t partial_dls, set<matrix_t> cms_set) {
    const short int n = partial_dls.size();
    assert(n > 0 and n < 11);
    matrix_t new_partial_dls(n, row_t(n,-1));
    matrix_t norm_partial_dls(n, row_t(n,-1));
    vector<matrix_t> dls_ascending_row;
    for (auto cms : cms_set) {
        for (short int i=0; i<n; i++) {
            for (short int j=0; j<n; j++) {
                const short int i2 = cms[i][j] / n;
                const short int j2 = cms[i][j] % n;
                new_partial_dls[i2][j2] = partial_dls[i][j];
            }
        }
        // Normalize by the first row:
        norm_partial_dls = normalize_first_row(new_partial_dls);
        dls_ascending_row.push_back(norm_partial_dls);
    }
    return dls_ascending_row;
}

set<matrix_t> find_main_class_repres_set(vector<matrix_t> dls_arr, set<matrix_t> cms_set) {
    set<matrix_t> main_class_repres_set;
    unsigned k=0;
    //unsigned max_main_class_repres_set_size = 0;
    for (auto &dls : dls_arr) {
        matrix_t main_class_repres = find_main_class_repres_all_cms(dls, cms_set);
        main_class_repres_set.insert(main_class_repres);
        /*
        if (main_class_repres_set.size() > max_main_class_repres_set_size) {
            max_main_class_repres_set_size = main_class_repres_set.size();
            cout << "main_class_repres_set size : " << main_class_repres_set.size() << endl;
        }
        */
        k++;
        if (k % 10000 == 0) cout << k << " out of " << dls_arr.size() << " processed" << endl;
    }
    return main_class_repres_set;
 }

int main(int argc, char *argv[]) {
    vector<string> argv_str;
    for (short int i=0; i<argc; i++) argv_str.push_back(argv[i]);
    if ((argc == 2) and (argv_str[1] == "-v")) {
        cout << prog << " of version " << version << endl;
        return 1;
    }
    if (argc < 2) {
        cout << "Usage : " << prog << " DLS-order [DLS-file] [CMS-file]" << endl;
        return 1;
    }
    short int n;
    istringstream(argv_str[1]) >> n;
    assert(n > 0 and n < 11);
    cout << prog << " of version " << version << " is running" << endl;
    cout << "DLS-order : " << n << endl;
    string dls_fname = "";
    string cms_fname = "";
    if (argc > 2) {
        dls_fname = argv_str[2];
        assert(dls_fname != "");
        cout << "DLS file name : " << dls_fname << endl;
    }
    if (argc > 3) {
        cms_fname = argv_str[3];
        assert(cms_fname != "");
        cout << "CMS file name : " << cms_fname << endl;
    }

    chrono::steady_clock::time_point start_point = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;

    matrix_t trivial_cms(n,row_t(n, 0));
    for (short int i=0; i<n; i++) {
        for (short int j=0; j<n; j++) {
            trivial_cms[i][j] = i*n + j;
        }
    }

    cout << "Trivial CMS :" << endl;
    print(trivial_cms);

    set<matrix_t> esodls_cms_set = reflect_rotate_symm(trivial_cms); 
    trivial_cms.clear();
    trivial_cms.shrink_to_fit();
    cout << "esodls_cms_set size : " << esodls_cms_set.size() << endl;

    cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;

    stringstream sstream_cms;
    sstream_cms << "cms_esodls_n" << n << "_new.txt";
    string out_cms_filename = sstream_cms.str();
    cout << "Writing ESODLS CMS to file " << out_cms_filename << endl;
    ofstream out_cms_file(out_cms_filename, ios_base::out);
    out_cms_file << "Writing " << esodls_cms_set.size() << " ESODLS CMS of order " << n << endl;
    unsigned k = 0;
    for (auto &cms : esodls_cms_set) {
        out_cms_file << "# " << k << endl;
        for (auto &row : cms) {
            for (short int j=0; j<row.size(); j++) {
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

    if (dls_fname == "" or cms_fname == "") {
        cout << "Stop since no DLS or no CMS are given" << endl;
        return 0;
    }

    vector<matrix_t> dls_arr = read_dls_string(dls_fname, n);
    cout << dls_arr.size() << " DLS were read" << endl;
    assert(dls_arr.size() > 0);

    set<matrix_t> cms_from_file = read_cms(cms_fname, n);
    cout << cms_from_file.size() << " CMS were read" << endl;
    assert(cms_from_file.size() > 0);

    set<matrix_t> diff_set;
    std::set_difference(
        cms_from_file.begin(), cms_from_file.end(),
        esodls_cms_set.begin(), esodls_cms_set.end(),
        std::inserter(diff_set, diff_set.end()));
    assert(cms_from_file.size() == esodls_cms_set.size() and diff_set.empty());
    cout << "CMS set is correct" << endl;
    cms_from_file.clear();
    diff_set.clear();

    cout << "Generating DLS main class" << endl;
    set<matrix_t> dls_main_class_repres_set = find_main_class_repres_set(dls_arr, esodls_cms_set);
    dls_arr.clear();
    dls_arr.shrink_to_fit();
    cout << "DLS main classes : " << dls_main_class_repres_set.size() << endl;
    stringstream sstream_main_class;
    sstream_main_class << "esodls_main_class_repres_n" << n << ".txt";
    string main_class_fname = sstream_main_class.str();
    cout << "Writing DLS main class representatives to file " << main_class_fname << endl;
    ofstream main_class_file(main_class_fname, ios_base::out);
    for (auto &x : dls_main_class_repres_set) {
        for (short int i=0; i<n; i++) {
            for (short int j=0; j<n; j++) {
                main_class_file << x[i][j];
                if (j != n - 1) main_class_file << " ";
            }
            main_class_file << endl;
        }
        main_class_file << endl;
    }
    main_class_file.close();

    cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;

    // Find all distinct ESODLS with ascending first row:
    cout << "Generaring DLS with ascending first row" << endl;
    set<vector<short int>> dls_ascending_row_set;
    k=0;
    for (auto &dls : dls_main_class_repres_set) {
        vector<matrix_t> dls_ascending_row = find_all_dls_with_ascending_row_from_main_class(dls, esodls_cms_set);
        for (auto &dls_asc_row : dls_ascending_row) {
            vector<short int> dls_usi_vec;
            short int usi = 0;
            unsigned digits_num = 0;
            for (auto &row : dls_asc_row) {
                for (auto &x : row) {
                    assert(x >= 0 and x < 11);
                    if (digits_num == 4) {
                        dls_usi_vec.push_back(usi);
                        digits_num = 0;
                        usi = 0;
                    }
                    usi += x*(short int)pow(10, digits_num);
                    digits_num++;
                }
            }
            dls_ascending_row_set.insert(dls_usi_vec);
        }
        k++;
        if (k % 1000 == 0) cout << k << " main classes out of " << dls_main_class_repres_set.size() << " processed" << endl;
    }
    cout << dls_ascending_row_set.size() << " DLS with ascending first row" << endl;
    dls_ascending_row_set.clear();

    cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;

    vector<matrix_t> x_based_basic_partial_dls = generate_x_based_basic_partial_dls(n);
    cout << x_based_basic_partial_dls.size() << " basic X-based fillings" << endl;
    cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;

    cout << "Finding main class representatives for X-based partial fillings" << endl;
    set<matrix_t> x_based_main_class_repres_set = find_main_class_repres_set(x_based_basic_partial_dls, esodls_cms_set);
    esodls_cms_set.clear();
    cout << "X-based partial filling main classes : " << x_based_main_class_repres_set.size() << endl;
    cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;
    stringstream sstream;
    sstream << "x_fillings_n" << n << "_new.txt";
    string x_fname = sstream.str();
    cout << "Writing X-based partial fillings to file " << x_fname << endl;
    ofstream ofile(x_fname, ios_base::out);
    k = 0;
    for (auto &x_fill : x_based_main_class_repres_set) {
        ofile << "X-based filling " << k << ":\n";
        for (short int i=0; i < n; i++) {
            for (short int j=0; j < n; j++) {
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

    cur_point = std::chrono::steady_clock::now();
    cout << "Elapsed " << std::chrono::duration_cast<std::chrono::seconds> (cur_point - start_point).count() << " sec" << std::endl;

    return 0;
}
