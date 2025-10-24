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
#include <thread>

#include "dlx_orth.h"

#include <omp.h>

#define time_point_t std::chrono::time_point<std::chrono::system_clock>

using namespace std;

string program = "dlx_mols";
string version = "0.2.1";

struct Record_orth_char_result {
	latinsquare_t square;
	unsigned transv;
	unsigned diag_transv;
};

struct Thread_result {
	vector<Record_orth_char_result> record_orth_char_results;
	unsigned max_orth_char;
	unsigned max_transv;
	unsigned num_max_transv;
	unsigned max_diag_transv;
	unsigned num_max_diag_transv;
};

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
vector<latinsquare_t> read_squares(const string filename, const unsigned n){
	vector<latinsquare_t> squares;
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
	unsigned long long num_processed_lines = 0;
	// The first mode - one row a line, elements divided by spaces:
	if (first_str.find(" ") != string::npos) {
		while (getline(in, s)) {
			if (s == "" and cur_square.size() == n){
				assert(DLX_orth::is_latinsquare(cur_square));
				squares.push_back(cur_square);
				cur_square.clear();
				num_processed_lines++;
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
			squares.push_back(cur_square);
			cur_square.clear();
			num_processed_lines++;
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
			squares.push_back(cur_square);
			cur_square.clear();
			num_processed_lines++;
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

// Find substring after a given prefix:
std::string str_after_prefix(const string str, const string prefix) {
    size_t pos = str.find( prefix );
    if ( pos != string::npos )
		return str.substr( pos + prefix.length() );
    return "";
}

void print_version() {
	cout << program << " of version " << version << endl;
}

void print_usage() {
	cout << "Usage : " << program << " DLS-file -cpunum=<int>" << endl;
	cout << "	-cpunum=<int> : (default = all cores) CPU cores\n";
}

int main(int argc, char *argv[]) 
{
	const unsigned n = 10;

	vector<string> str_argv;
	for (int i=0; i < argc; ++i) str_argv.push_back(argv[i]);
	assert(str_argv.size() == argc);

	if (argc == 2 and str_argv[1] == "-v") {
		print_version();
		exit(EXIT_SUCCESS);
	}
	if (argc == 2 and str_argv[1] == "-h") {
		print_usage();
		exit(EXIT_SUCCESS);
	}
  	if (argc < 3) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	string filename = argv[1];
	unsigned cpu_num = 0;
	if (str_argv.size() > 2) {
		string s = str_after_prefix(str_argv[2], "-cpunum=");
		if (s != "") istringstream(s) >> cpu_num;
	}
	cout << "cpu_num	: " << cpu_num << endl;
	// If no CPU num is given, use all threads:
	if (!cpu_num) {
		cpu_num = thread::hardware_concurrency();
		cout << cpu_num << " CPU cores are detected" << endl;
	}
	omp_set_num_threads(cpu_num);
	cout << cpu_num << " threads are used" << endl;

	const time_point_t program_start = std::chrono::system_clock::now();

	vector<latinsquare_t> squares = read_squares(filename, n);
	cout << squares.size() << " Latin squares were read" << endl;
	cout << current_time(program_start) << endl;
	unsigned num_diagls = 0;
	for (auto &square : squares) if (DLX_orth::is_diag_latinsquare(square)) num_diagls++;
	cout << "Of them " << num_diagls << " diagonal Latin squares" << endl;
	cout << current_time(program_start) << endl;
	cout << "Start processing squares" << endl;

	set<latinsquare_t> max_orth_char_squares;

	unsigned num_squares_74orth = 0;
	unsigned max_orth_char = 0;
	unsigned max_diag_transv = 0;
	unsigned max_transv = 0;
	vector<Thread_result> threads_results;
	threads_results.resize(cpu_num);
	for (auto &t_r : threads_results) {
		t_r.max_orth_char = 0;
		t_r.max_transv = 0;
		t_r.num_max_transv = 0;
		t_r.max_diag_transv = 0;
		t_r.num_max_diag_transv = 0;
	}

	unsigned k = 0;
	unsigned report_per_task = 10000;
	if (cpu_num >= 10) report_per_task = 100000;
	if (cpu_num >= 100) report_per_task = 1000000;
	// There are plenty of simple tasks, so static distribution sounds here:
	//#pragma omp parallel for schedule(static, 1000)
	#pragma omp parallel for schedule(static, 1000)
	for (auto &square : squares) {
		if ((k % report_per_task == 0) && (k > 0)) cout << k << " squares processed" << endl;
		LS_result ls_res = DLX_orth::find_transversals_and_orth_mates(square);
		// Increase the loop-counter by all threads, but not at the same time:
		#pragma omp critical
		{
			k++;
		}
		unsigned transv = ls_res.transv;
		unsigned diag_transv = ls_res.diag_transv;
		assert(transv > 0);
		assert(diag_transv > 0);
		int thread_id = omp_get_thread_num();
		// Process the number of transv:
		if (threads_results[thread_id].max_transv == 0 or transv > threads_results[thread_id].max_transv) {
			// The thread-wise record is updated:
			threads_results[thread_id].max_transv = transv;
			// Set the counter of squares with the record number of transversals to 0:
			threads_results[thread_id].num_max_transv = 0;
			// Check if a global record is updated as well (one thread a time):
			#pragma omp critical
			if (max_transv == 0 or transv > max_transv) {
				max_transv = transv;
				cout << "New maximum number of tranvsersals : " << max_transv << endl;
			}
		}
		// Increase the number of squares with the record number of transversals:
		if (transv == threads_results[thread_id].max_transv) {
			threads_results[thread_id].num_max_transv++;
		}
		// Process the number of diagonal transversals:
		if (threads_results[thread_id].max_diag_transv == 0 or diag_transv > threads_results[thread_id].max_diag_transv) {
			// The thread-wise record is updated:
			threads_results[thread_id].max_diag_transv = diag_transv;
			// Set the counter to 0:
			threads_results[thread_id].num_max_diag_transv = 0;
			// Check if a global record is updated as well:
			#pragma omp critical
			if (max_diag_transv == 0 or diag_transv > max_diag_transv) {
				max_diag_transv = diag_transv;
				cout << "New maximum number of diagonal tranvsersals : " << diag_transv << endl;
			}
		}
		// Increase the number of squares with the record transversals number:
		if (diag_transv == threads_results[thread_id].max_diag_transv) {
			threads_results[thread_id].num_max_diag_transv++;
		}
		
		// For all pairs of DLS which are orthogonal to the current square and form a triple:
		if (ls_res.orth_mates.size() < 2) continue;
		for (unsigned j = 0; j < ls_res.orth_mates.size(); j++) {
			assert(DLX_orth::is_diag_latinsquare(ls_res.orth_mates[j]));
		}
		for (unsigned j = 0; j < ls_res.orth_mates.size() - 1; j++) {
			for (unsigned j2 = j+1; j2 < ls_res.orth_mates.size(); j2++) {
				unsigned orth_char = calc_orth_char(ls_res.orth_mates[j], ls_res.orth_mates[j2]);
				if (orth_char == 74) {
					#pragma omp critical
					num_squares_74orth++;
				}
				// Check if the thread-wise record is updated:
				if (threads_results[thread_id].max_orth_char == 0 or orth_char > threads_results[thread_id].max_orth_char) {
					threads_results[thread_id].max_orth_char = orth_char;
					threads_results[thread_id].record_orth_char_results.clear();
					// If the global maximum orth char is updated:
					#pragma omp critical
					if (max_orth_char == 0 or orth_char > max_orth_char) {
						max_orth_char = orth_char;
						cout << "Updated max_orth_char : " << max_orth_char << endl;
						cout << current_time(program_start) << endl;
					}
				}
				// Save squares with the thread-wise record orth char:
				if (orth_char == max_orth_char) {
					Record_orth_char_result res;
					res.square = square;
					res.transv = transv;
					res.diag_transv = diag_transv;
					threads_results[thread_id].record_orth_char_results.push_back(res);
				} 
				assert(orth_char <= max_orth_char);
			}
		}
	}

	// Collect record squares:
	cout << "max_transv : " << max_transv << endl;
	cout << "max_diag_transv : " << max_diag_transv << endl;
	vector<Record_orth_char_result> record_squares;
	for (auto &x : threads_results) {
		if (x.max_orth_char == max_orth_char) {
			for (auto &y : x.record_orth_char_results) {
				record_squares.push_back(y);
			}
		}
	}
	assert(record_squares.size() == num_squares_74orth);
	cout << record_squares.size() << " squares with maximum orthogonal char " << max_orth_char << endl;
	// Write the record squares to a file:
	stringstream sstream;
	sstream << "orth_char=" << max_orth_char << "_squares_order=10";
	string max_orth_char_squares_fname = sstream.str();
	cout << "Writing " << record_squares.size() << 
		" squares to file " << max_orth_char_squares_fname << endl;
	ofstream ofile(max_orth_char_squares_fname, ios_base::out);
	ofile << "transv diag_transv square" << endl;
	for (auto &res : record_squares) {
		ofile << res.transv << " " << res.diag_transv << " ";
		for (unsigned i = 0; i < res.square.size(); i++) {
			assert(not res.square[i].empty());
			for (unsigned j = 0; j < res.square[i].size(); j++) {
				ofile << res.square[i][j];
			}
		}
		ofile << endl;
	}
	ofile.close();

	cout << "End" << endl;
	cout << current_time(program_start) << endl;

	return 0;
}
