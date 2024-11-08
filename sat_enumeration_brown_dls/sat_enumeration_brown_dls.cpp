// Created on: 6 Nov 2024
// Author: Oleg Zaikin
// E-mail: oleg.zaikin@icc.ru
//
// For a given n, enumerate all Brown diagonal Latin squares of order n
// via an AllSAT solver.
//
// Example:
//   ./sat_enumeration_brown_dls 8 ./cms_n8 -cpunum=10
//=============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <string>
#include <sstream>
#include <chrono>
#include <thread>
#include <cassert>
#include <cmath>

#include <omp.h>

#define clause_t vector<int>
#define cell_t vector<unsigned>
#define cover_t vector<vector<int>>
#define row_t vector<short int>
#define matrix_t vector<row_t>
#define time_point_t std::chrono::time_point<std::chrono::system_clock>

using namespace std;

string program = "sat_enumeration_brown_dls";
string version = "0.1.3";

struct SatEncDls {
    vector<vector<cell_t>> X;
    unsigned vars_num;
    vector<clause_t> clauses;
};

string clause_to_str(const clause_t &cla);
vector<clause_t> at_most_one_clauses(vector<unsigned> &vars);
vector<clause_t> exactly_one_clauses(vector<unsigned> &vars);
vector<clause_t> latin_square_clauses(const vector<vector<cell_t>> &X, const bool &is_diag);
SatEncDls diag_latin_square(const unsigned &n, const bool &is_first_row_ascending);
vector<clause_t> equality_clauses(const unsigned &x, const unsigned &y);
vector<clause_t> vertical_symmetry_clauses(const vector<vector<cell_t>> &X,
                                           const unsigned &n, const cover_t &cover);
vector<clause_t> horizontal_symmetry_clauses(const vector<vector<cell_t>> &X,
                                             const unsigned &n);
int conseq_multip(const int &lower_bound, const int &upper_bound);
vector<vector<int>> combinations(const int &n, const int &k);
vector<cover_t> gen_covers(const unsigned &n);
string replace(const string &str, const string &orig_substr, const string &repl_substr);
string generate_cnf_brown_dls_vertic_sym(const unsigned &n, const SatEncDls &satencdls,
                                         const cover_t &cover, const unsigned &cover_index);
string generate_cnf_brown_dls_horiz_sym(const unsigned &n, const SatEncDls &satencdls,
                                        const cover_t &cover, const unsigned &cover_index);
string ls_from_sat(const vector<int> &sat_assignment, const unsigned &n);
vector<string> parse_latin_squares_from_sat_assign(const string &solver_out_fname,
                                                   const unsigned &n);
vector<string> find_all_dls_sat_solver(const string &cnf_fname, const unsigned &n,
                                       const unsigned &cover_index);
bool is_digits(const string &str);
set<matrix_t> read_cms(const string &filename, const unsigned &n);
string min_main_class_repres(const string &dls, const unsigned &n, const set<matrix_t> &cms_set);
string normalize(const string &ls, const unsigned &n);
void print(const cover_t &cover);
void print(const string &ls);
string str_after_prefix(const string &str, const string &prefix);
string current_time(const time_point_t &program_start);

int main(int argc, char *argv[])
{
    vector<string> argv_str;
    for (unsigned i=0; i<argc; i++) argv_str.push_back(argv[i]);
	if ((argc == 2) and (argv_str[1] == "-v")) {
        cout << program << " of version " << version << endl;
        return 1;
    }

    if (argc < 4) {
        cout << "Usage: " << program << " DLS-order CMS-file -cpunum=<int>" << endl;
        return 1;
    }

    const time_point_t program_start = std::chrono::system_clock::now();

    unsigned n = 0;
    istringstream(argv_str[1]) >> n;
    assert(n > 1 and n < 11 and n % 2 == 0);
    cout << "Running " << program << " of version " << version << endl;
    cout << "n : " << n << endl;
    string cms_fname = "";
    set<matrix_t> cms_set;
    
    cms_fname = argv_str[2];
    cout << "CMS file name : " << cms_fname << endl;
    cms_set = read_cms(cms_fname, n);
    assert(cms_set.size() > 0);
    cout << cms_set.size() << " CMS were read" << endl;

    string s = str_after_prefix(argv_str[3], "-cpunum=");
    unsigned cpunum = 0;
	if (s != "") std::istringstream(s) >> cpunum;

    const unsigned nthreads = cpunum > 0 ? cpunum : std::thread::hardware_concurrency();
	std::cout << "threads       : " << nthreads << std::endl;
	omp_set_num_threads(nthreads);

    // Make a CNF that encodes a diagonal Latin square of order n:
    SatEncDls satencdls = diag_latin_square(n, true);
    cout << satencdls.vars_num << " variables encode DLS" << endl;
    cout << satencdls.clauses.size() << " clauses encode DLS" << endl;

    string cnf_dls_fname = "dls_n" + to_string(n) + "_first_row_ascending.cnf";
    //cout << "Writing a CNF for DLS to the file " << cnf_dls_fname << endl;
    ofstream ofile(cnf_dls_fname, ios_base::out);
    ofile << "p cnf " + to_string(satencdls.vars_num) + " " + to_string(satencdls.clauses.size()) + "\n";
    for (auto &cla : satencdls.clauses) ofile << clause_to_str(cla) + '\n';
    ofile.close();

    vector<cover_t> covers = gen_covers(n);

    // For each cover, generate a horizontal-symmetry and a vertical-symmetry CNF:
    set<string> main_class_repres_set;
    unsigned processed = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (unsigned cover_index=0; cover_index < covers.size(); cover_index++) {
        cout << "Before generating horiz CNFs " << current_time(program_start) << endl;
        string horiz_sym_cnf_fname = generate_cnf_brown_dls_horiz_sym(n, satencdls, covers[cover_index], cover_index);
        vector<string> ls_str_arr_horiz = find_all_dls_sat_solver(horiz_sym_cnf_fname, n, cover_index);
        cout << "After parsing all horiz DLS from sat assignments " << current_time(program_start) << endl;
        if (cms_set.size() > 0) {
            for (auto &ls : ls_str_arr_horiz) {
                string min_repres = min_main_class_repres(ls, n, cms_set);
                assert(min_repres.size() == n*n);
                main_class_repres_set.insert(min_repres);
            }
        }
        cout << "After finding all horiz main class representatives " << current_time(program_start) << endl;
        //
        string vertic_sym_cnf_fname = generate_cnf_brown_dls_vertic_sym(n, satencdls, covers[cover_index], cover_index);
        vector<string> ls_str_arr_vertic = find_all_dls_sat_solver(vertic_sym_cnf_fname, n, cover_index);
        if (cms_set.size() > 0) {
            for (auto &ls : ls_str_arr_vertic) {
                string min_repres = min_main_class_repres(ls, n, cms_set);
                assert(min_repres.size() == n*n);
                main_class_repres_set.insert(min_repres);
            }
        }
        if (ls_str_arr_horiz.size() > 0 or ls_str_arr_vertic.size() > 0) {
            string brown_dls_cover_fname = "brown_dls_n" + to_string(n) + "_cover" + to_string(cover_index);
            ofstream ofile(brown_dls_cover_fname, ios_base::out);
            for (auto &ls : ls_str_arr_horiz) ofile << ls << endl;
            for (auto &ls : ls_str_arr_vertic) ofile << ls << endl;
            ofile.close();
        }
        processed++;
        cout << "processed " << processed << " covers out of " << covers.size() << endl;
        cout << main_class_repres_set.size() << " main class representatives" << endl;
        cout << current_time(program_start) << endl;
    }

    cout << main_class_repres_set.size() << " main class representatives" << endl;
    string main_class_fname = "brown_dls_n" + to_string(n) + "_main_class_repres";
    ofstream main_class_file(main_class_fname, ios_base::out);
    for (auto &dls : main_class_repres_set) {
        for (unsigned i=0; i<n; i++) {
            for (unsigned j=0; j<n; j++) {
                main_class_file << dls[i*n + j];
                if (j != n-1) main_class_file << " ";
            }
            main_class_file << endl;
        }
        main_class_file << endl;
    }
    main_class_file.close();

    return 0;
}

string current_time(const time_point_t &program_start) {
    stringstream sstream;
    const time_point_t program_end = std::chrono::system_clock::now();
    sstream << "Elapsed : "
    << std::chrono::duration_cast<std::chrono::seconds>(program_end - program_start).count()
    << " seconds";
    return sstream.str();
}

string clause_to_str(const clause_t &cla) {
    assert(not cla.empty());
    string str = "";
    for (auto &lit : cla) str += to_string(lit) + " ";
    return str + "0";
}

// Clauses for the AtMostOne constraint:
vector<clause_t> at_most_one_clauses(vector<unsigned> &vars) {
	assert(vars.size() > 1);
	vector<clause_t> res_clauses;
	for (unsigned i=0; i<vars.size()-1; i++) {
		for (unsigned j=i+1; j<vars.size(); j++) {
            assert(i != j);
            clause_t cla = {-(int)vars[i], -(int)vars[j]};
            res_clauses.push_back(cla);
        }
    }
	return res_clauses;
}

// Clauses for the ExactlyOne constraint for a set of variables:
vector<clause_t> exactly_one_clauses(vector<unsigned> &vars) {
	assert(vars.size() > 1);
    vector<clause_t> res_clauses;
    // At most one constraint: 
	res_clauses = at_most_one_clauses(vars);
    // At Least One constraint:
    clause_t alo_clause;
    for (auto &v : vars) alo_clause.push_back((int)v);
	res_clauses.push_back(alo_clause); 
	return res_clauses;
}

vector<clause_t> latin_square_clauses(const vector<vector<cell_t>> &X, const bool &is_diag) {
    const unsigned n = X.size();
    assert(n > 0 and n < 11);
    vector<clause_t> dls_clauses;
    vector<clause_t> eo_clauses;
	// Constraints on rows, columns, and values are obligatory:
	for (unsigned i=0; i<n; i++) {
		for (unsigned j=0; j<n; j++) {
			// Each square's cell contains exactly one value 0..n-1:
            vector<unsigned> cell_variables;
            for (unsigned k=0; k<n; k++) cell_variables.push_back(X[i][j][k]);
			eo_clauses = exactly_one_clauses(cell_variables);
            for (auto &cla : eo_clauses) dls_clauses.push_back(cla);
			// Each value occurs exactly once in each row:
            vector<unsigned> row_variables;
            for (unsigned k=0; k<n; k++) row_variables.push_back(X[i][k][j]);
			eo_clauses = exactly_one_clauses(row_variables);
            for (auto &cla : eo_clauses) dls_clauses.push_back(cla);
			// Each value occurs exactly once in each column:
            vector<unsigned> col_variables;
            for (unsigned k=0; k<n; k++) col_variables.push_back(X[k][i][j]);
			eo_clauses = exactly_one_clauses(col_variables);
            for (auto &cla : eo_clauses) dls_clauses.push_back(cla);
        }
    }
	// Constraints on main diagonal and antidiagonal are optional:
	if (is_diag) {
		// Main diagonal:
		for (unsigned i=0; i<n; i++) {
            vector<unsigned> main_diag_variables;
            for (unsigned k=0; k<n; k++) main_diag_variables.push_back(X[k][k][i]);
			eo_clauses = exactly_one_clauses(main_diag_variables);
            for (auto &cla : eo_clauses) dls_clauses.push_back(cla);
        }
		// Main antidiagonal:
        for (unsigned i=0; i<n; i++) {
            vector<unsigned> main_antidiag_variables;
            for (unsigned k=0; k<n; k++) main_antidiag_variables.push_back(X[n-1-k][k][i]);
			eo_clauses = exactly_one_clauses(main_antidiag_variables);
            for (auto &cla : eo_clauses) dls_clauses.push_back(cla);
        }
    }
    return dls_clauses;
}

SatEncDls diag_latin_square(const unsigned &n, const bool &is_first_row_ascending) {
    assert(n > 0 and n < 11);
    // n^3 variables in the CNF encode a Latin square:
    vector<vector<cell_t>> X(n, vector<cell_t>(n, cell_t(n, 0)));
    // Variables for Latin squares X:
    unsigned vars_num = 0;
    // The first variable is 1:
    for (unsigned i=0; i<n; i++) {
        for (unsigned j=0; j<n; j++) {
            for (unsigned k=0; k<n; k++) X[i][j][k] = ++vars_num;
        }
    }
    // The first row in ascending order:
    vector<clause_t> clauses;
    if (is_first_row_ascending) {
        for (unsigned j=0; j<n; j++) {
            clause_t cla = {(int)X[0][j][j]};
            clauses.push_back(cla);
        }
    }
    // clauses for diagonal Latin square:
    vector<clause_t> dls_clauses = latin_square_clauses(X, true);
    assert(dls_clauses.size() > 0);
    for (auto &cla : dls_clauses) clauses.push_back(cla);
    assert(dls_clauses.size() > 0);
    SatEncDls satencdls;
    satencdls.X = X;
    assert(vars_num == n*n*n);
    satencdls.vars_num = vars_num;
    satencdls.clauses = clauses;
    return satencdls;
}

// Clauses for equality of two variables:
vector<clause_t> equality_clauses(const unsigned &x, const unsigned &y) {
    assert(x > 0 and y > 0 and x != y);
    vector<clause_t> clauses;
    clause_t cla1 = {(int)x, -(int)y};
    clause_t cla2 = {-(int)x, (int)y};
    clauses.push_back(cla1);
    clauses.push_back(cla2);
    return clauses;
}

vector<clause_t> vertical_symmetry_clauses(const vector<vector<cell_t>> &X, const unsigned &n, const cover_t &cover) {
    assert(n > 0 and n < 11);
    assert(n % 2 == 0 and X.size() == n);
    vector<clause_t> res_clauses;
    unsigned half_n = (unsigned)(n/2);
    assert(cover.size() == half_n);
    for (unsigned i=0; i<half_n; i++) { // for each row 0 .. (n/2)-1
        for (unsigned j=0; j<n; j++) { // for each column
            for (auto &inverse_columns_indices : cover) {
                assert(inverse_columns_indices.size() == 2);
                vector<clause_t> eq_clauses = equality_clauses(X[i][j][inverse_columns_indices[0]], X[n-1-i][j][inverse_columns_indices[1]]);
                for (auto &cla : eq_clauses) res_clauses.push_back(cla);
            }
        }
    }
    return res_clauses;
}

vector<clause_t> horizontal_symmetry_clauses(const vector<vector<cell_t>> &X, const unsigned &n) {
    assert(n > 0 and n < 11);
    assert(n % 2 == 0 and X.size() == n);
    vector<clause_t> res_clauses;
    unsigned half_n = (unsigned)(n/2);
    for (unsigned i=0; i<n; i++) { // for each row
        for (unsigned j=0; j<half_n; j++) { // for each column 0 .. (n/2)-1
            for (unsigned k=0; k<n; k++) { // for each cell's value
                // x[i][j] == n-1-x[i][n-1-j]
                vector<clause_t> eq_clauses = equality_clauses(X[i][j][k], X[i][n-1-j][n-1-k]);
                for (auto &cla : eq_clauses) res_clauses.push_back(cla);
            }
        }
    }
    return res_clauses;
}

// Multiplication lower_bound * (lower_bound + 1) * ... * upper_bound
int conseq_multip(const int &lower_bound, const int &upper_bound) {
    assert(lower_bound > 0 and upper_bound > 0 and lower_bound <= upper_bound);
	int final_val = 1;
	for (int i = lower_bound; i <= upper_bound; i++) final_val *= i;
	return final_val;
}

// Generate all k-combinations of n elements:
vector<vector<int>> combinations(const int &n, const int &k) {
    assert(n > 0 and k > 0 and k <= n);
    vector<vector<int>> combs;
	int val;
	vector<int> index_arr;
	index_arr.resize(k);
	unsigned comb_index = 0;
	combs.resize(conseq_multip(n - k + 1, n) / conseq_multip(1, k));
	for ( unsigned i = 0; i < combs.size(); i++ )
		combs[i].resize(k);
	for ( unsigned i = 0; i < index_arr.size(); i++ )
		index_arr[i] = i; // start indexes
	for ( ;; ) {
		combs[comb_index++] = index_arr;
		if ( index_arr[k-1] != n-1 ) { // increase last value if we can
			index_arr[k-1]++;
			continue;
		}
		val = k-1;
		while ( ( val >= 0 ) and ( index_arr[val] == (n-k) + val ) )
			val--; // val - index of rightest value that is not on final position

		if ( val >= 0 ) { // if some values on final positions but not all
			index_arr[val]++;
			for ( int i = val+1; i < k; i++ )
				index_arr[i] = index_arr[i-1] + 1; // set initial state to all values right to changed one
		} // if val < 0 then final state - all '1' in tail
		else
			break; // all values are on final positions
	}
    return combs;
}

void print(const string &ls) {
    const unsigned ls_len = ls.size();
    assert(ls_len > 0);
    const unsigned n = (unsigned)sqrt(ls_len);
    assert(n > 0);
    for (unsigned i=0; i<n; i++) {
        for (unsigned j=0; j<n; j++) {
            assert(i*n + j < ls_len);
            cout << ls[i*n + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print(const cover_t &cover) {
    assert(cover.size() > 0);
    for (auto &pair : cover) {
        assert(pair.size() == 2);
        cout << "(" << pair[0] << "," << pair[1] << ") ";
    }
    cout << endl;
}

vector<cover_t> gen_covers(const unsigned &n) {
    unsigned num_row_pairs = (unsigned)(n*(n-1)/2);
    cout << "num_row_pairs : " << num_row_pairs << endl;
    vector<vector<int>> all_pairs = combinations(n, 2);
    /*
    for (auto &comb : pair_combinations) {
        for (auto &x : comb) cout << x << " ";
        cout << endl;
    }
    */
    assert(all_pairs.size() == num_row_pairs);
    int half_n = (int)(n/2);
    vector<vector<int>> pairs_combinations_indices = combinations(all_pairs.size(), half_n);
    vector<vector<vector<int>>> pairs_combinations;
    for (auto &indices : pairs_combinations_indices) {
        vector<vector<int>> cur_pairs;
        for (auto &i : indices) cur_pairs.push_back(all_pairs[i]);
        pairs_combinations.push_back(cur_pairs);
    }
    cout << pairs_combinations.size() << " pairs combinations" << endl;
    // Here a cover is a set of n/2 pairs which contain all n rows:
    vector<cover_t> covers;
    unsigned comb_num = 0;
    for (auto &comb : pairs_combinations) {
        comb_num += 1;
        set<int> row_indices_set;
        for (auto &pair : comb) {
            row_indices_set.insert(pair[0]);
            row_indices_set.insert(pair[1]);
        }
        assert(row_indices_set.size() <= n);
        if (row_indices_set.size() != n) continue;
        else covers.push_back(comb);
    }
    assert(covers.size() > 0);
    cout << covers.size() << " covers of " << half_n << " pairs out of " << comb_num << " are possible" << endl;
    cout << "The first covers are :" << endl;
    print(covers[0]);
    if (covers.size() > 1) {
        print(covers[1]);
    }
    if (covers.size() > 2) {
        cout << "The last cover is :" << endl;
        print(covers[covers.size()-1]);
    }
    return covers;
}

string replace(const string &str, const string &orig_substr, const string &repl_substr) {
    assert(not str.empty());
    size_t pos = str.find(orig_substr);
    if (pos == string::npos) return str;
    return str.substr(0, pos) + repl_substr + str.substr(pos+orig_substr.size(), str.size()-pos-orig_substr.size());
}

string generate_cnf_brown_dls_vertic_sym(const unsigned &n,
                                         const SatEncDls &satencdls,
                                         const cover_t &cover,
                                         const unsigned &cover_index) {
    assert(n > 0 and n % 2 == 0);
    assert(satencdls.vars_num > 0);
    assert(satencdls.X.size() > 0);
    assert(satencdls.clauses.size() > 0);
    assert(cover.size() == (unsigned)(n/2));
    // Add constraints for horizontal symmetry:
    vector<clause_t> vertic_sym_clauses = vertical_symmetry_clauses(satencdls.X, n, cover);
    //cout << vertic_sym_clauses.size() << " clauses encode vertical symmetry" << endl;
    assert(vertic_sym_clauses.size() > 0);
    // Inverse-columns for a given cover:
    vector<clause_t> all_eq_clauses;
    for (auto &inverse_columns_indices : cover) {
        assert(inverse_columns_indices.size() == 2);
        for (unsigned i=0; i<n; i++) {
            for (unsigned k=0; k<n; k++) {
                vector<clause_t> eq_clauses = equality_clauses(satencdls.X[i][inverse_columns_indices[0]][k],
                                                               satencdls.X[n-1-i][inverse_columns_indices[1]][k]);
                for (auto &cla : eq_clauses) all_eq_clauses.push_back(cla);
            }
        }
    }
     // Write CNF to file:
    string cnf_name = "dls_n" + to_string(n) + "_first_row_ascending_vertic_sym_inverse_rows_" + to_string(cover_index) + ".cnf";
    //cout << "Writing CNF for DLS with vertical symmetry and inverse columns to file " << cnf_name << endl;
    const unsigned cla_num = satencdls.clauses.size() + vertic_sym_clauses.size() + all_eq_clauses.size();
    ofstream ofile(cnf_name, ios_base::out);
    ofile << "p cnf " + to_string(satencdls.vars_num) + " " + to_string(cla_num) + "\n";
    for (auto &cla : satencdls.clauses) ofile << clause_to_str(cla) + "\n";
    for (auto &cla : vertic_sym_clauses) ofile << clause_to_str(cla) + "\n";
    for (auto &cla : all_eq_clauses) ofile << clause_to_str(cla) + "\n";
    ofile.close();
    return cnf_name;
}

string generate_cnf_brown_dls_horiz_sym(const unsigned &n,
                                        const SatEncDls &satencdls,
                                        const cover_t &cover,
                                        const unsigned &cover_index) {
    assert(n > 0 and n % 2 == 0);
    assert(satencdls.vars_num > 0);
    assert(satencdls.X.size() > 0);
    assert(satencdls.clauses.size() > 0);
    assert(cover.size() == (unsigned)(n/2));
    // Add constraints for horizontal symmetry:
    vector<clause_t> horiz_sym_clauses = horizontal_symmetry_clauses(satencdls.X, n);
    assert(horiz_sym_clauses.size() > 0);
    //cout << horiz_sym_clauses.size() << " clauses encode horizontal symmetry" << endl;
    // Add row-inverse constraints:
    vector<clause_t> all_eq_clauses;
    for (auto &inverse_rows_indices : cover) {
        assert(inverse_rows_indices.size() == 2);
        for (unsigned j=0; j<n; j++) {
            for (unsigned k=0; k<n; k++) {
                vector<clause_t> eq_clauses = equality_clauses(satencdls.X[inverse_rows_indices[0]][j][k],
                                                               satencdls.X[inverse_rows_indices[1]][n-1-j][k]);
                for (auto &cla : eq_clauses) all_eq_clauses.push_back(cla);
            }
        }
    }
    // Write CNF to file:
    string cnf_name = "dls_n" + to_string(n) + "_first_row_ascending_horiz_sym_inverse_rows_" + to_string(cover_index) + ".cnf";
    //cout << "Writing CNF for DLS with horizontal symmetry and inverse rows to file " << cnf_name << endl;
    const unsigned cla_num = satencdls.clauses.size() + horiz_sym_clauses.size() + all_eq_clauses.size();
    ofstream ofile(cnf_name, ios_base::out);
    ofile << "p cnf " + to_string(satencdls.vars_num) + " " + to_string(cla_num) + "\n";
    for (auto &cla : satencdls.clauses) ofile << clause_to_str(cla) + "\n";
    for (auto &cla : horiz_sym_clauses) ofile << clause_to_str(cla) + "\n";
    for (auto &cla : all_eq_clauses) ofile << clause_to_str(cla) + "\n";
    ofile.close();
    return cnf_name;
}

// Given a SAT assignment, form the LS:
string ls_from_sat(const vector<int> &sat_assignment, const unsigned &n) {
	assert(sat_assignment.size() > 0);
	assert(n > 0);
	assert(sat_assignment.size() == n*n*n);
	string ls_str = "";
	vector<int> cell_literals;
	for (auto &lit : sat_assignment) {
		// Literals of a current cell are collected:
		cell_literals.push_back(lit);
		int val = -1;
		unsigned neg_lit_num = 0;
		if (cell_literals.size() == n) {
			for (unsigned j=0; j < n; ++j) {
				if (cell_literals[j] > 0) {
					assert(val == -1);
					val = (int)j;
				}
				else neg_lit_num += 1;
			}
			assert(neg_lit_num == n - 1);
			assert(val >= 0 and (unsigned)val <= n - 1);
			ls_str += to_string(val);
			cell_literals.clear();
		}
	}
	assert(ls_str.size() == n * n);
	return ls_str;
}

// Parse satysfying assignments from a solver's output and convert them to Latin squares:
vector<string> parse_latin_squares_from_sat_assign(const string &solver_out_fname, const unsigned &n) {
    assert(n > 0 and n < 11);
    assert(not solver_out_fname.empty());
	string line;
    ifstream f(solver_out_fname, ios_base::in);
    
	bool is_started_sol = false;
	vector<int> cur_assignment;
    vector<string> ls_str_arr;

	// Process all lines:
	while (getline(f, line)) {
		if (line.size() < 2) continue;
		if ( (line.find("s SATISFIABLE") != string::npos) or
		     (line.find("c Answer: ") != string::npos) )
		{
            is_started_sol = true;
            cur_assignment.clear();
            continue;
		}
		if (is_started_sol) {
            // If clasp solver, no solution after after 's SATISFIABLE':
            if (line[0] != 'v' or line[1] != ' ') continue;
            // Cut 'v' at the beginning:
            line.erase(line.begin());
            // Read literals
            stringstream sstream(line);
            vector<int> literals;
            int ival;
            while (sstream >> ival) literals.push_back(ival);
            cur_assignment.insert(cur_assignment.end(), literals.begin(), literals.end());
            // If the last line of a solution is found, stop reading:
            if (cur_assignment[cur_assignment.size() - 1] == 0) {
                is_started_sol = false;
                // Cut '0' at the end:
                cur_assignment.pop_back();
                // Check that there are N^3 variables:
                //assert(cur_assignment.size() == 3*pow(ls_order,3));
                // Add the assignment to the list:
                string ls_str = ls_from_sat(cur_assignment, n);
                assert(ls_str.size() == n*n);
                ls_str_arr.push_back(ls_str);
            }
		}
	}

    f.close();
    cout << ls_str_arr.size() << " DLS were parsed" << endl;
    return ls_str_arr;
}

vector<string> find_all_dls_sat_solver(const string &cnf_fname, const unsigned &n, const unsigned &cover_index) {
    assert(n > 0 and n < 11 and n % 2 == 0);
    // Find all satisfying assignments:
    string solver_out_fname = "out_" + replace(cnf_fname, ".cnf", "") + "_cover_" + to_string(cover_index);
    string sys_str = "clasp --models 0 --enum-mode=bt --configuration=crafty ./" + cnf_fname + " > " + solver_out_fname;
    int systemRet = system(sys_str.c_str());
    assert(systemRet >= 0);
    // Parse Latin squares from sat assignments:
    vector<string> ls_str_arr = parse_latin_squares_from_sat_assign(solver_out_fname, n);
    sys_str = "rm " + cnf_fname;
    systemRet = system(sys_str.c_str());
    assert(systemRet == 0);
    sys_str = "rm " + solver_out_fname;
    systemRet = system(sys_str.c_str());
    assert(systemRet == 0);
    return ls_str_arr;
}

bool is_digits(const string &str)
{
    return all_of(str.begin(), str.end(), ::isdigit);
}

// Read CMSs from a given file:
set<matrix_t> read_cms(const string &filename, const unsigned &n) {
    assert(n > 0 and n < 11);
    assert(not filename.empty());
	set<matrix_t> cms_set;
	ifstream in;
	in.open(filename);
	assert(in.is_open());
	string s;
	matrix_t cms;
	while (getline(in, s)) {
		if (s == "" and cms.size() > 0){
            assert(cms.size() == n);
			cms_set.insert(cms);
			cms.clear();
		}
		else {
			vector<short int> tmp;
            stringstream sstream(s);
            string word;
            while (sstream >> word) {
                if (not is_digits(word)) {
                    tmp.clear();
                    break;
                }
                short int si = -1;
                istringstream(word) >> si;
				tmp.push_back(si);
            }
            if (not tmp.empty()) {
                assert(tmp.size() == n);
                cms.push_back(tmp);
            }
		}
	}
	in.close();
	return cms_set;
}

// Normalize LS, i.e. make its main diagonal 0, 1, ..., n-1
string normalize(const string &ls, const unsigned &n) {
    assert(n > 0 and n < 11);
    string norm_ls(n*n, '-');
    string norm_perm(n, '-');
    // Element in permutation is i if its index is i-th element in main diag:
    for (unsigned i = 0; i < n; i++) {
        // Here ls[i*n + i] - '0' is char->int convertion to get index:
        norm_perm[ls[i*n + i] - '0'] = (char)i + '0';
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            norm_ls[i*n + j] = norm_perm[ls[i*n + j] - '0'];
        }
    }
    for (unsigned i = 0; i < n; i++) {
        assert(norm_ls[i*n + i] == (char)i + '0');
    }
    return norm_ls;
}

// Given a DLS as a string and a set of CMS, find the minimal main class representative:
string min_main_class_repres(const string &dls, const unsigned &n, const set<matrix_t> &cms_set) {
    assert(n > 0 and n < 11);
    assert(dls.size() == n*n);
    assert(cms_set.size() > 0);
    string new_dls(n*n, '-');
    string norm_dls = "";
    string min_repres = "";
    for (auto cms : cms_set) {
        assert(cms.size() == n);
        for (unsigned i=0; i<n; i++) {
            assert(cms[i].size() == n);
            for (unsigned j=0; j<n; j++) {
                const unsigned i2 = cms[i][j] / n;
                const unsigned j2 = cms[i][j] % n;
                assert(i2*n + j2 < new_dls.size());
                assert(new_dls.size() == dls.size());
                new_dls[i2*n + j2] = dls[i*n + j];
            }
        }
        // Normalize:
        norm_dls = normalize(new_dls, n);
        if (min_repres == "" or norm_dls < min_repres) {
            min_repres = norm_dls;
        }
    }
    assert(min_repres.size() == n*n);
    return min_repres;
}

string str_after_prefix(const string &str, const string &prefix) {
    int found = str.find( prefix );
    if ( found != -1 )
	return str.substr( found + prefix.length() );
    return "";
}
