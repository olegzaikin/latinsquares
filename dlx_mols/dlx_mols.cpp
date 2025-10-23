// Created on: 15 Oct 2024
// Author: Oleg Zaikin
// E-mail: oleg.zaikin@icc.ru
//
// For a given file with Latin squares, find all their orthogonal mates
// by the DLX algorithm.
//
// Example:
//   ./dlx_mols ./squares
//==========================================================================

#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include <cassert>
#include <map>

#include "dlx_orth.h"

#define time_point_t std::chrono::time_point<std::chrono::system_clock>

using namespace std;

string program = "dlx_mols";
string version = "0.1.14";

int strtoi(string s) {
	assert(not s.empty());
	int x = atoi(s.c_str());
	return x;
}

string inttostr(int number) {
	stringstream sstream;
	sstream << number;
	return sstream.str();
}

string current_time(const time_point_t &program_start) {
    stringstream sstream;
    const time_point_t program_end = std::chrono::system_clock::now();
    sstream << "Elapsed : "
    << std::chrono::duration_cast<std::chrono::seconds>(program_end - program_start).count()
    << " seconds";
    return sstream.str();
}

// Print a given Latin square:
void print_square(const latinsquare_t square) {
	assert(not square.empty());
	for (unsigned i = 0; i < square.size(); i++) {
		assert(not square[i].empty());
		for (unsigned j = 0; j < square[i].size(); j++) {
			cout << square[i][j] << " ";
		}
		cout << endl;
	}
}

// Read Latin squares from a given file:
set<latinsquare_t> read_squares(const string filename, const unsigned n){
	set<latinsquare_t> squares_set;
	latinsquare_t cur_square;
	ifstream in;
	string s;
	in.open(filename);
	assert(in.is_open());
	// Read the first line to choose a reading mode:
	string first_str;
	getline(in, first_str);
	in.close();
	// Close the file and open again - now for reading squares:
	in.open(filename);
	unsigned long long processed_lines_num = 0;
	// The first mode - one row a line, elements divided by spaces:
	if (first_str.find(" ") != string::npos) {
		while (getline(in, s)) {
			if (s == "" and cur_square.size() == n){
				assert(DLX_orth::is_latinsquare(cur_square));
				squares_set.insert(cur_square);
				cur_square.clear();
				processed_lines_num++;
			}
			else if (s == "") continue;
			else {
				row_t row;
				for (unsigned j = 0; j < n; j++) {
					string t_s = s.substr(2 * j, 1);
					row.push_back(strtoi(t_s));
				}
				cur_square.push_back(row);
			}
		}
		if (cur_square.size() == n) {
			squares_set.insert(cur_square);
			cur_square.clear();
			processed_lines_num++;
		}
	}
	// The second mode - one square a line without spaces:
	else {
		while (getline(in, s)) {
			if (s == "") continue;
			assert(s.size() == n*n);
			for (unsigned i=0; i<n; i++) {
				row_t row;
				for (unsigned j=0; j<n; j++) {
					assert(i*n + j < s.size());
					// One char to string:
					string t_s(1, s[i*n + j]);
					// String to int:
					row.push_back(strtoi(t_s));
				}
				assert(row.size() == n);
				cur_square.push_back(row);
			}
			assert(cur_square.size() == n);
			squares_set.insert(cur_square);
			cur_square.clear();
			processed_lines_num++;
		}
	}
	in.close();
	cout << processed_lines_num << " squares were read" << endl; 
	return squares_set;
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
	cout << "Running " << program << " of version " << version << endl;
	if (argc < 2) {
        cout << "Usage : " << program << " DLS-file" << endl;
        return 1;
    }

	string filename = argv[1];
	const unsigned n = 10;

	const time_point_t program_start = std::chrono::system_clock::now();

	set<latinsquare_t> squares_set = read_squares(filename, n);
	cout << squares_set.size() << " Latin unique squares were read" << endl;
	cout << current_time(program_start) << endl;
	unsigned diagls_num = 0;
	for (auto &square : squares_set) if (DLX_orth::is_diag_latinsquare(square)) diagls_num++;
	cout << "Of them " << diagls_num << " diagonal Latin squares" << endl;
	cout << current_time(program_start) << endl;
	cout << "Start processing squares" << endl;

	set<latinsquare_t> max_orth_char_squares;
	// A key is the number of diagonal transversals, while a value is the number of squares:
	map<unsigned, unsigned> diag_transversals_num_sizes;
	map<unsigned, unsigned> transversals_num_sizes;

	unsigned max_orth_char = 0;
	unsigned max_diag_transv_num = 0;
	unsigned max_transv_num = 0;
	unsigned k = 0;
	for (auto &square : squares_set) {
		if ((k % 10000 == 0) && (k > 0)) cout << k << " squares processed" << endl;
		LS_result ls_res = DLX_orth::find_transversals_and_orth_mates(square);
		unsigned cur_diag_transv_num = ls_res.diag_transversals_num;
		assert(cur_diag_transv_num > 0);
		if (diag_transversals_num_sizes.find(cur_diag_transv_num) == diag_transversals_num_sizes.end()) {
			diag_transversals_num_sizes.insert(pair<unsigned,unsigned>(cur_diag_transv_num,0));
			if (max_diag_transv_num == 0 or cur_diag_transv_num > max_diag_transv_num) {
				max_diag_transv_num = cur_diag_transv_num;
				cout << "New maximum number of diagonal tranvsersals : " << cur_diag_transv_num << endl;
			}
		}
		unsigned cur_transv_num = ls_res.transversals_num;
		assert(cur_transv_num > 0);
		if (transversals_num_sizes.find(cur_transv_num) == transversals_num_sizes.end()) {
			transversals_num_sizes.insert(pair<unsigned,unsigned>(cur_transv_num,0));
			if (max_transv_num == 0 or cur_transv_num > max_transv_num) {
				max_transv_num = cur_transv_num;
				cout << "New maximum number of tranvsersals : " << cur_transv_num << endl;
			}
		}
		// Increase the transversals-sizes counter:
		diag_transversals_num_sizes[cur_diag_transv_num]++;
		transversals_num_sizes[cur_transv_num]++;
		k++;
		// cout << orth_mates.size() << endl;
		// For all pairs of DLS which are orthogonal to the current square and form a triple:
		if (ls_res.orth_mates.size() < 2) continue;
		for (unsigned j = 0; j < ls_res.orth_mates.size(); j++) {
			assert(DLX_orth::is_diag_latinsquare(ls_res.orth_mates[j]));
		}
		for (unsigned j = 0; j < ls_res.orth_mates.size() - 1; j++) {
			for (unsigned j2 = j+1; j2 < ls_res.orth_mates.size(); j2++) {
				unsigned orth_char = calc_orth_char(ls_res.orth_mates[j], ls_res.orth_mates[j2]);
				// If max orth char is updated:
				if (max_orth_char == 0 or orth_char > max_orth_char) {
					max_orth_char = orth_char;
					cout << "Updated max_orth_char : " << max_orth_char << endl;
					cout << ls_res.orth_mates.size() << " orthogonal mates" << endl;
					// New record, so clear all squares from the previous record:
					max_orth_char_squares.clear();
				}
				assert(orth_char <= max_orth_char);
				// If one more squares with the current max orth char is found:
				if (orth_char != 0 and orth_char == max_orth_char) {
					cout << "***" << endl;
					print_square(square);
					max_orth_char_squares.insert(square);
					cout << max_orth_char_squares.size() << " squares with orth char " << max_orth_char << endl;
					cout << cur_diag_transv_num << " diagonal transversals" << endl;
					cout << cur_transv_num << " transversals" << endl;
					cout << current_time(program_start) << endl;
					cout << endl;
				}
			}
		}
	}

	stringstream sstream;
	sstream << "orth_char=" << max_orth_char << "_squares_order=10";
	string max_orth_char_squares_fname = sstream.str();
	cout << "Writing " << max_orth_char_squares.size() << 
		" squares to file " << max_orth_char_squares_fname << endl;

	ofstream ofile(max_orth_char_squares_fname, ios_base::out);
	for (auto &square : max_orth_char_squares) {
		for (unsigned i = 0; i < square.size(); i++) {
			assert(not square[i].empty());
			for (unsigned j = 0; j < square[i].size(); j++) {
				ofile << square[i][j] << " ";
			}
			ofile << endl;
		}
		ofile << endl;
	}
	ofile.close();

	string out_diag_transversals_sizes_fname = "transversals-diag_num_sizes_order=10";
	cout << "Writing sizes of diagonal transverals to file " << out_diag_transversals_sizes_fname << endl;
	ofstream ofile2(out_diag_transversals_sizes_fname, ios_base::out);
	ofile2 << "transversals-diag : squares" << endl;
	for (auto pair : diag_transversals_num_sizes) {
		ofile2 << pair.first << " : " << pair.second << endl;
	}
	ofile2.close();

	string out_transversals_sizes_fname = "transversals_num_sizes_order=10";
	cout << "Writing sizes of transverals to file " << out_transversals_sizes_fname << endl;
	ofstream ofile3(out_transversals_sizes_fname, ios_base::out);
	ofile3 << "transversals : squares" << endl;
	for (auto pair : transversals_num_sizes) {
		ofile3 << pair.first << " : " << pair.second << endl;
	}
	ofile3.close();

	cout << "End" << endl;
	cout << current_time(program_start) << endl;

	return 0;
}
