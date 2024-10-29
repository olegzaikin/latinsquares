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
#include <cassert>

#include "dlx_orth.h"

using namespace std;

string program = "dlx_mols";
string version = "0.1.3";

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

// Print a given Latin square:
void print_square(const latinsquare_t square) {
	assert(not square.empty());
	cout << endl;
	for (unsigned i = 0; i < square.size(); i++) {
		assert(not square[i].empty());
		for (unsigned j = 0; j < square[i].size(); j++) {
			cout << square[i][j] << " ";
		}
		cout << endl;
	}
}

// Read Latin squares from a given file:
vector<latinsquare_t> read_squares(const string filename, const unsigned n){
	vector<latinsquare_t> squares;
	ifstream in;
	in.open(filename);
	assert(in.is_open());
	string s;
	latinsquare_t cur_square;
	while (getline(in, s)) {
		if (s == "" and cur_square.size() > 0){
			assert(DLX_orth::is_latinsquare(cur_square));
			squares.push_back(cur_square);
			cur_square.clear();
		}
		else if (s == "") continue;
		else {
			vector<int> tmp;
			for (unsigned j = 0; j < n; j++) {
				string t_s = s.substr(2 * j, 1);
				tmp.push_back(strtoi(t_s));
			}
			cur_square.push_back(tmp);
		}
	}
	in.close();
	return squares;
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

	vector<latinsquare_t> squares = read_squares(filename, n);
	cout << squares.size() << " Latin squares were read" << endl;
	unsigned diagls_num = 0;
	for (auto &square : squares) if (DLX_orth::is_diag_latinsquare(square)) diagls_num++;
	cout << "Of them " << diagls_num << " diagonal Latin squares" << endl;
	
	unsigned max_orth_char = 0;
	for (int i = 0; i < squares.size(); i++) {
		if ((i % 10000 == 0) && (i > 0)) cout << i << " squares processed" << endl;
		vector<latinsquare_t> orth_mates = DLX_orth::find_all_orth_mates(squares[i]);
		// cout << orth_mates.size() << endl;
		// For all pairs of DLS which are orthogonal to the current square and form a triple:
		if (orth_mates.size() < 2) continue;
		for (unsigned j = 0; j < orth_mates.size() - 1; j++) {
			for (unsigned j2 = j+1; j2 < orth_mates.size(); j2++) {
				unsigned orth_char = calc_orth_char(orth_mates[j], orth_mates[j2]);
				//if (orth_char == max_orth_char) cout << orth_char << endl;
				if (max_orth_char == 0 or orth_char > max_orth_char) {
					max_orth_char = orth_char;
					cout << "Updated max_orth_char : " << max_orth_char << endl;
					print_square(squares[i]);
				}
			}
		}
	}

	return 0;
}
