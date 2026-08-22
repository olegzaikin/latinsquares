// Created on: 22 Aug 2026
// Author: Oleg Zaikin
// E-mail: oleg.zaikin@icc.ru
//
// For a given n, ESODLS CMS, and triples of DLS,
// find all main classes for the first DLSs from the triples
//
// Example:
//   ./check_main_class_dls_for_triples 10 ./cms_esodls_n10_new.txt ./orth_char=74_squares_order=10_brown_dls_n10_horiz_nodupl
//=============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <sstream>
#include <cassert>
#include <cmath>

#define row_t vector<short int>
#define matrix_t vector<row_t>

using namespace std;

string program = "check_main_class_dls_for_triples";
string version = "0.0.1";

set<matrix_t> read_cms(const string &filename, const unsigned &n);
vector<string> read_first_dlss_from_triples(const string triples_fname, const unsigned n);
string min_main_class_repres(const string &dls, const unsigned &n, const set<matrix_t> &cms_set);
string normalize(const string &ls, const unsigned &n);
bool is_digits(const string &str);

int main(int argc, char *argv[])
{
    vector<string> argv_str;
    for (unsigned i=0; i<argc; i++) argv_str.push_back(argv[i]);
	if ((argc == 2) and (argv_str[1] == "-v")) {
        cout << program << " of version " << version << endl;
        return 1;
    }

    if (argc < 4) {
        cerr << "Usage: " << program << " DLS-order CMS-file TriplesDLS" << endl;
        return 1;
    }

    unsigned n = 0;
    istringstream(argv_str[1]) >> n;
    if ( ((n < 2) or (n > 10)) or (n % 2 != 0) ) {
        cerr << "Error. n must be >= 2, <= 10, and even." << endl;
        exit(1);
    }
    cout << "Running " << program << " of version " << version << endl;
    cout << "n : " << n << endl;
    string cms_fname = "";
    set<matrix_t> cms_set;
    
    cms_fname = argv_str[2];
    cout << "CMS file name : " << cms_fname << endl;
    cms_set = read_cms(cms_fname, n);
    if ( cms_set.empty() ) {
        cerr << "Error. CMS set is empty." << endl;
        exit(1);
    }
    cout << cms_set.size() << " CMS were read" << endl;

    string triples_fname = argv_str[3];
    cout << "Triples file name : " << triples_fname << endl;
    vector<string> first_dlss_str = read_first_dlss_from_triples(triples_fname, n);
    if ( first_dlss_str.empty() ) {
        cerr << "Error. No first DLSs from triples." << endl;
        exit(1);
    }
    cout << first_dlss_str.size() << " first DLSs from triples were read" << endl;

    set<string> min_repres_set;
    for (unsigned i=0; i < first_dlss_str.size(); i++) {
        if (i > 0 and i % 1000 == 0) cout << i << " DLSs are processed" << endl;
        string min_repres = min_main_class_repres(first_dlss_str[i], n, cms_set);
        min_repres_set.insert(min_repres);
    }

    cout << min_repres_set.size() << " main classes of diagonal Latin squares" << endl;
    cout << *min_repres_set.begin() << endl;

    return 0;
}

vector<string> read_first_dlss_from_triples(const string triples_fname, const unsigned n) {
    vector<string> first_dlss_str;
    assert(n > 0 and n < 11);
    assert(not triples_fname.empty());
	set<matrix_t> cms_set;
	ifstream in;
	in.open(triples_fname);
	assert(in.is_open());
	string s;
    // Skipe the first line:
    getline(in, s);
	while (getline(in, s)) {
        if (s == "") continue;
        vector<short int> tmp;
        stringstream sstream(s);
        string word1, word2, word3;
        sstream >> word1;
        sstream >> word2;
        sstream >> word3;
        assert(word3.size() == n*n);
        first_dlss_str.push_back(word3);
        //cout << word3;
        /*
        matrix_t dls;
        row_t row;
        short int si;        
        for (auto &x : word3) {
            si = static_cast<short int>(x);
            row.push_back(si);
            if (row.size() == n) {
                assert(row.size() == n);
                dls.push_back(row);
                row.clear();
            }
        }
        assert(dls.size() == n);
        first_dlss.push_back(dls);
        */
    }

    return first_dlss_str;
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


// Normalize LS by the first row, i.e. make it 0, 1, ..., n-1
string normalize(const string &ls, const unsigned &n) {
    assert(n > 0 and n < 11);
    string norm_ls(n*n, '-');
    string norm_perm(n, '-');
    // Element in permutation is i if its index is i-th element in main diag:
    for (unsigned i = 0; i < n; i++) {
        // Here ls[i*n + i] - '0' is char->int convertion to get index:
        norm_perm[ls[i] - '0'] = (char)i + '0';
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            norm_ls[i*n + j] = norm_perm[ls[i*n + j] - '0'];
        }
    }
    for (unsigned i = 0; i < n; i++) {
        assert(norm_ls[i] == (char)i + '0');
    }
    return norm_ls;
}
