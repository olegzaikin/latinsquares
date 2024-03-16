// Created on: 15 March 2024
// Author: Oleg Zaikin
// E-mail: zaikin.icc@gmail.com
//
// XXX
//
// Usage : XXX
//
// Example:
//     XXX
//=============================================================================

#include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <set>

#pragma warning(disable : 4996)

using namespace std;

string program = "solve_sat_cms_mpi";
string version = "0.0.1";

#define cms_t vector<vector<unsigned>>
#define partial_ls_t vector<vector<string>>

vector<cms_t> cms_vec;
vector<partial_ls_t> partial_ls_vec;

enum wu_status{ NOT_STARTED = -1, IN_PROGRESS = 0, PROCESSED = 1 };
enum wu_result{ UNKNOWN = 0, UNSAT = 1, SAT = 2, INTERR = 3 };

const int REPORT_EVERY_SEC = 100;

struct workunit {
	int id;
	int cms_index;
	int partial_ls_index;
	wu_status status;
	wu_result result;
	double processing_time;
	workunit() : id(-1), cms_index(-1), partial_ls_index(-1), 
		status(NOT_STARTED), result(UNKNOWN), processing_time(-1) {};
};

struct CNF {
	long long int var_num;
	long long int clause_num;
	vector<string> clauses;
	CNF() : var_num(0), clause_num(0), clauses() {}
	CNF(string cnf_name) : var_num(0), clause_num(0), clauses() {
		read(cnf_name);
	}
	void read(string cnf_name) {
		ifstream cnf_file(cnf_name, ios_base::in);
		string str;
		while (getline(cnf_file, str)) {
			if (str.size() == 0 or str[0] == 'p' or str[0] == 'c')
				continue;
			clauses.push_back(str);
			clause_num++;
			stringstream sstream;
			sstream << str;
			long long int ival;
			while (sstream >> ival)	var_num = max(llabs(ival), var_num);
		}
		cnf_file.close();
	}
	void print() {
		cout << "var_num : " << var_num << endl;
		cout << "clause_num : " << clause_num << endl;
	};
};

// MPI functions:
void control_process(const int corecount, const string cms_file_name, 
	const string partial_ls_file_name);
void computing_process(const int rank, const string solver_file_name,
                       const string cnf_file_name,
											 const string cms_file_name, 
											 const string partial_ls_file_name);
void send_wu(workunit &wu, const int wu_index, const int computing_process_id);
// Output functions:
void write_info_out_file(const string control_process_ofile_name,
                      vector<workunit> &wu_vec, const double start_time);
void write_processing_info(vector<workunit> &wu_vec);
// Latin-square-based functions:
vector<workunit> generate_workunits(const string cms_file_name, 
																		const string partial_ls_file_name);
vector<cms_t> read_cms(const string cms_file_name);
vector<partial_ls_t> read_partial_ls(
	const string partial_ls_file_name,
  const unsigned ls_order);
unsigned cnf_var_num(const unsigned ls_order, const unsigned ls_index,
                     const unsigned row_index, const unsigned col_index,
										 const unsigned cell_val);
wu_result find_all_mols(string solver_name, const CNF cnf, workunit &wu);
vector<vector<int>> parse_sat_assignments( const string output_str,
	const bool is_file_output, const unsigned ls_order);
string first_ls_from_sat(const vector<int> sat_assignment,
	const unsigned ls_order);
void enumerate_by_allsat_solver(string solver_name,
														    const string cur_cnf_name,
														    const string local_out_file_name);
// Other functions:
string str_after_prefix(const string str, const string prefix);
string exec(const string cmd_str);

void print_usage() {
	cout << "Usage : ./" << program
						<< " solver CNF CMS partial-LS\n";
	cout << " solver     	  : AllSAT solver\n";
  cout << " CNF        	  : file with CNF that encodes 2MOLS\n";
  cout << " CMS        	  : file with list of CMSs\n";
  cout << " partial-LS 	  : file with list of partial LSs\n";
}

void print_version() {
	cout << program << " of version " << version << endl;
}

int main(int argc, char *argv[])
{
	int rank = 0, corecount = 0;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	vector<string> str_argv;
	for (unsigned i=0; i<argc; ++i)
		str_argv.push_back(argv[i]);
	assert((unsigned)argc == str_argv.size());

	if (rank == 0) {
		cout << "Running " << program << " of version " << version << endl;
		if ( (argc == 2) && 
		     ((str_argv[1] == "-v") || (str_argv[1] == "--version")) ) {
		    return 1;
		}
		cout << "corecount " << corecount << endl;
	}

	if ( (argc < 5) && (rank == 0) ) {
		print_usage();
		return 1;
	}

  assert(corecount > 1);
	assert(rank >= 0);

	string solver_file_name     = str_argv[1];
	string cnf_file_name	      = str_argv[2];
	string cms_file_name        = str_argv[3];
	string partial_ls_file_name = str_argv[4];

	// Control process:
	if (rank == 0) {
		cout << "solver_name          : " << solver_file_name   << endl;
		cout << "cnf_file_name        : " << cnf_file_name << endl;
		cout << "cms_file_name        : " << cms_file_name << endl;
		cout << "partial_ls_file_name : " << partial_ls_file_name << endl;
		control_process(corecount, cms_file_name, partial_ls_file_name);
	}
	else { 
		// Computing process:
		computing_process(rank, solver_file_name, cnf_file_name, cms_file_name, 
		  partial_ls_file_name);
	}

	return 0;
}


vector<cms_t> read_cms(const string cms_file_name) {
	vector<cms_t> cms_arr;
	ifstream cms_file(cms_file_name, ios_base::in);
	assert(cms_file.is_open());
	string str;
	unsigned ls_order = 0;
	cms_t cms;
	while (getline(cms_file, str)) {
		if ( (str.size() < 2) or (str[0] == '#') ) continue;
		if ( (str.find("Loading") != string::npos) or
			   (str.find("ESODLS") != string::npos)
				) continue;
		vector<unsigned> row;
		stringstream sstream;
		sstream << str;
		unsigned uval;
		while (sstream >> uval) row.push_back(uval);
		if (!ls_order) ls_order = row.size();
		assert(ls_order);
		cms.push_back(row);
		if (cms.size() == ls_order) {
			cms_arr.push_back(cms);
			cms.clear();
		}
	}
	cms_file.close();

	return cms_arr;
}

vector<partial_ls_t> read_partial_ls( const string partial_ls_file_name,
  const unsigned ls_order)
{
	vector<partial_ls_t> partial_ls_arr;
	assert(ls_order);

	ifstream partial_ls_file(partial_ls_file_name, ios_base::in);
	assert(partial_ls_file.is_open());
	string str;
	partial_ls_t partial_ls;
	while (getline(partial_ls_file, str)) {
		if ( (str.size() < 2) or 
		     (str.find("filling") != string::npos) 
			 )
				continue;
		vector<string> row;
		stringstream sstream;
		sstream << str;
		string word;
		while (sstream >> word) row.push_back(word);
		assert(row.size() == ls_order);
		partial_ls.push_back(row);
		if (partial_ls.size() == ls_order) {
			assert(partial_ls.size() == ls_order);
			partial_ls_arr.push_back(partial_ls);
			partial_ls.clear();
		}
	}
	partial_ls_file.close();

	return partial_ls_arr;
}

// Generate workunits:
vector<workunit> generate_workunits(const string cms_file_name, 
																		const string partial_ls_file_name) 
{
	// Read all CMSs from file:
	cms_vec = read_cms(cms_file_name);
	assert(cms_vec.size() > 0);
	cout << cms_vec.size() << " CMSs were read" << endl;
	cout << "The first CMS is :" << endl;
	for (auto x : cms_vec[0]) {
		for (auto y : x) {
			cout << y << " ";
		}
		cout << endl;
	}
	unsigned ls_order = cms_vec[0].size();
	cout << endl;
	cout << "LS order : " << ls_order << " was detected in CMS" << endl;
	cout << endl;

	partial_ls_vec = read_partial_ls(partial_ls_file_name, ls_order);
	cout << partial_ls_vec.size() << " partial Latin squares were read" << endl;
	cout << "The first one is" << endl;
	for (auto x : partial_ls_vec[0]) {
		for (auto y : x) {
			cout << y << " ";
		}
		cout << endl;
	}

	vector<workunit> wu_vec;
	unsigned wu_id = 0;
	for (unsigned i=0; i<cms_vec.size(); ++i) {
		for (unsigned j=0; j<partial_ls_vec.size(); ++j) {
			workunit wu;
			wu.id = wu_id;
			wu.cms_index = i;
			wu.partial_ls_index = j;
			wu.status = NOT_STARTED;
			wu.result = UNKNOWN;
			wu.processing_time = -1;
			wu_vec.push_back(wu);
			wu_id++;
		}
	}
	assert(wu_vec.size() > 0);
	assert(wu_vec.size() == cms_vec.size() * partial_ls_vec.size());

	return wu_vec;
}

// A process that generates and manages tasks for computing processes:
void control_process(const int corecount, const string cms_file_name,
	const string partial_ls_file_name)
{
	double start_time = MPI_Wtime();
	
	vector<workunit> wu_vec = generate_workunits(cms_file_name, partial_ls_file_name);
	cout << wu_vec.size() << " workunits were generated" << endl;
	
	// Erase progress file:
	string control_process_ofile_name = "!total_progress";
	ofstream control_process_ofile(control_process_ofile_name, ios_base::out);
	control_process_ofile.close();
	
	// Send the first portion of workunits:
	int wu_index = 0;
	for (int i = 0; i < corecount - 1; i++) {
		assert(wu_index >= 0 and wu_index < wu_vec.size());
		send_wu(wu_vec[wu_index], wu_index, i + 1);
		wu_index++;
		// If wus num is lower than core count:
		if (wu_index == wu_vec.size()) break;
	}

	// Receive results and send back new workunits:
	unsigned total_processed_wus = 0;
	int processed_wu_index;
	wu_result res;
	double time;
	int stop_mes = -1;
	MPI_Status status, current_status;
	double sum_runtime = 0.0;
	double result_writing_time = -1;
	while (total_processed_wus < wu_vec.size()) {
		// Receive a result:
		MPI_Recv(&processed_wu_index, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		current_status = status;
		MPI_Recv(&res,  1, MPI_INT, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&time, 1, MPI_DOUBLE, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
		assert((res >= UNSAT) && (res <= INTERR));
		assert(wu_vec[processed_wu_index].status == IN_PROGRESS);
		wu_vec[processed_wu_index].status = PROCESSED;
		wu_vec[processed_wu_index].result = res;
		wu_vec[processed_wu_index].processing_time = time;
		sum_runtime += time;
		total_processed_wus++;
	
		// Send back a new workunit if one exists:
		if (wu_index < wu_vec.size()) {
			assert(wu_index >= 0 and wu_index < wu_vec.size());
			send_wu(wu_vec[wu_index], wu_index, current_status.MPI_SOURCE);
			wu_index++;
		}
		// Send the stop-message:
		else {
			cout << "sending stop message to computing process " << current_status.MPI_SOURCE << endl;
			cout << "total_processed_wus : " << total_processed_wus << endl;
			MPI_Send(&stop_mes, 1, MPI_INT, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}

		if ( (result_writing_time < 0) || 
		     (MPI_Wtime() - result_writing_time > REPORT_EVERY_SEC) )
		{
			write_info_out_file(control_process_ofile_name, wu_vec, start_time);
			write_processing_info(wu_vec);
			result_writing_time = MPI_Wtime();
		}
	}

	write_info_out_file(control_process_ofile_name, wu_vec, start_time);
	write_processing_info(wu_vec);
	cout << "control process has finished" << endl;
	
	MPI_Finalize();
	//MPI_Abort(MPI_COMM_WORLD, 0);
}

// Send a task from the control process to a computing process:
void send_wu(workunit &wu, const int wu_index, const int computing_process_id)
{
	MPI_Send(&wu_index, 1, MPI_INT, computing_process_id, 0, MPI_COMM_WORLD);
	assert(wu.status == NOT_STARTED);
	wu.status = IN_PROGRESS;
}

// Write output data to a file:
void write_info_out_file(const string control_process_ofile_name,
                         vector<workunit> &wu_vec, const double start_time)
{
	assert(wu_vec.size() > 0);
	assert(start_time > 0.0);
	double min_solving_time_unsat = 1e+308;
	double max_solving_time_unsat = -1;
	double avg_solving_time_unsat = -1;
	double sum_time_unsat = 0.0;
	double min_solving_time_sat = 1e+308;
	double max_solving_time_sat = -1;
	double avg_solving_time_sat = -1;
	double sum_time_sat = 0.0;
	unsigned sat_sol = 0;
	unsigned unsat_sol = 0;
	unsigned processed = 0;

	for (auto cur_wu : wu_vec) {
		assert(cur_wu.result != INTERR);
		if (cur_wu.status != PROCESSED) continue;
		if (cur_wu.result == UNSAT) {
			unsat_sol++;
			max_solving_time_unsat = max(cur_wu.processing_time, max_solving_time_unsat);
			min_solving_time_unsat = min(cur_wu.processing_time, min_solving_time_unsat);
			sum_time_unsat += cur_wu.processing_time;
		}
		else if (cur_wu.result == SAT) {
			sat_sol++;
			max_solving_time_sat = max(cur_wu.processing_time, max_solving_time_sat);
			min_solving_time_sat = min(cur_wu.processing_time, min_solving_time_sat);
			sum_time_sat += cur_wu.processing_time;
		}
		processed++;
	}
	if (sum_time_unsat > 0)
		avg_solving_time_unsat = sum_time_unsat / unsat_sol;
	if (sum_time_sat > 0)
		avg_solving_time_sat = sum_time_sat / sat_sol;
	
	double percent_val;
	ofstream control_process_ofile(control_process_ofile_name, ios_base::app);
	control_process_ofile << endl << "***" << endl;
	control_process_ofile << "elapsed time  : "
		<< MPI_Wtime() - start_time << " sec." << endl;
	control_process_ofile << "total WUs     : " << wu_vec.size() << endl;
	percent_val = double(processed * 100) / (double)wu_vec.size();
	control_process_ofile << "processed wus : " << processed
		<< ", i.e. " << percent_val << " %" << endl;
	control_process_ofile << "* sat   : " << sat_sol << endl;
	control_process_ofile << "mintime : " << min_solving_time_sat << endl;
	control_process_ofile << "maxtime : " << max_solving_time_sat << endl;
	control_process_ofile << "avgtime : " << avg_solving_time_sat << endl;
	control_process_ofile << "* unsat : " << unsat_sol << endl;
	control_process_ofile << "mintime : " << min_solving_time_unsat << endl;
	control_process_ofile << "maxtime : " << max_solving_time_unsat << endl;
	control_process_ofile << "avgtime : " << avg_solving_time_unsat << endl;
	control_process_ofile << endl;
	control_process_ofile.close();
}

// Write info about all tasks:
void write_processing_info(vector<workunit> &wu_vec)
{
	ofstream ofile("!processing_info");
	ofile << "id result time" << endl;
	for (auto &cur_wu : wu_vec)
		ofile << cur_wu.id << " " << cur_wu.result << " " << cur_wu.processing_time << endl;
	ofile.close();
	ofile.clear();
}

void computing_process(const int rank, const string solver_file_name,
                       const string cnf_file_name,
											 const string cms_file_name, 
											 const string partial_ls_file_name)
{
	// Form workunits once:
	vector<workunit> wu_vec = generate_workunits(cms_file_name, partial_ls_file_name);
	// Read CNF once:
	CNF cnf(cnf_file_name);

	assert(wu_vec.size() > 0);
	assert(cms_vec.size() > 0);
	assert(partial_ls_vec.size() > 0);
	
	MPI_Status status;
	int wu_index = -1;
	int wu_id = -1;
	for (;;) {
		MPI_Recv( &wu_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
		//cout << "received wu_index " << wu_index << endl;
		if (wu_index == -1) {// stop message
			cout << "computing prosess " << rank << " got the stop message" << endl;
			break;
		}

		double start_time = MPI_Wtime();
		wu_result res = find_all_mols(solver_file_name, cnf, wu_vec[wu_index]);
		double solving_time = MPI_Wtime() - start_time;

		// send calculated result to the control process
		//cout << "sending wu_index " << wu_index << endl;
		//cout << "sending res " << res << endl;
		MPI_Send( &wu_index,     1, MPI_INT,    0, 0, MPI_COMM_WORLD);
		MPI_Send( &res,          1, MPI_INT,    0, 0, MPI_COMM_WORLD);
		MPI_Send( &solving_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	} // for (;;)

	MPI_Finalize();
}

// Find substring after a given prefix:
string str_after_prefix(const string str, const string prefix) {
    size_t pos = str.find( prefix );
    if ( pos != string::npos )
			return str.substr( pos + prefix.length() );
    return "";
}

// Given information on a LS's cell, return its variable in a CNF:
unsigned cnf_var_num(const unsigned ls_order, const unsigned ls_index,
                     const unsigned row_index, const unsigned col_index,
										 const unsigned cell_val)
{
	return ls_index*pow(ls_order,3) + row_index*pow(ls_order,2)
		+ col_index*ls_order + cell_val + 1;
}

// Find all MOLS in a CNF with added CMS- and partial-LS-constraints:
wu_result find_all_mols(string solver_name, const CNF cnf, workunit &wu)
{
	unsigned ls_order = partial_ls_vec[wu.partial_ls_index].size();
	assert(ls_order > 0);

	// CMS constraints:
	vector<string> clauses_cms;
  for (unsigned i=0; i < ls_order; ++i) {
		for (unsigned j=0; j < ls_order; ++j) {
			vector<unsigned> first_ls_cell_vars;
			// CNF variables for values in cell (i,j) of LS 1
			for (unsigned k=0; k < ls_order; ++k) {
				first_ls_cell_vars.push_back(cnf_var_num(ls_order, 0, i, j, k));
			}
			assert(first_ls_cell_vars.size() == ls_order);
			// cell (i, j) in LS 1 == cell (i2, j2) in LS 2:
			unsigned i2 = cms_vec[wu.cms_index][i][j] / ls_order;
      unsigned j2 = cms_vec[wu.cms_index][i][j] % ls_order;
			vector<unsigned> second_ls_cell_vars;
			// CNF variables for values in cell (i2,j2) of LS 2
			for (unsigned k=0; k < ls_order; ++k) {
				second_ls_cell_vars.push_back(cnf_var_num(ls_order, 1, i2, j2, k));
			}
			assert(second_ls_cell_vars.size() == ls_order);
			// Equality-clauses for each pair of values:
			for (unsigned k=0; k < ls_order; ++k) {
				string str_val1 = to_string(first_ls_cell_vars[k]);
				string str_val2 = to_string(second_ls_cell_vars[k]);
				clauses_cms.push_back(      str_val1 + " -" + str_val2 + " 0");
				clauses_cms.push_back("-" + str_val1 + " "  + str_val2 + " 0");
			}
		}
	}
	assert(clauses_cms.size() == pow(ls_order,3)*2);

	// Partial-LS constraints for the first LS:
	vector<int> partial_ls_literals;
	vector<string> clauses_partial_ls;
	for (unsigned i=0; i < ls_order; ++i) {
		for (unsigned j=0; j < ls_order; ++j) {
			string val_str = partial_ls_vec[wu.partial_ls_index][i][j];
			if (val_str == "-") continue;
			unsigned val_uint;
			istringstream(val_str) >> val_uint;
			unsigned var = cnf_var_num(ls_order, 0, i, j, val_uint);
			clauses_partial_ls.push_back(to_string(var) + " 0");
			partial_ls_literals.push_back(var);
		}
	}

	string cur_cnf_name = "./cms" + to_string(wu.cms_index) + "_" +
		"partial-ls" + to_string(wu.partial_ls_index) + ".cnf";

	// Write all needed clauses to a file:
	ofstream ofile(cur_cnf_name, ios_base::out);
	unsigned clauses_num = cnf.clauses.size() + clauses_partial_ls.size() +
		clauses_cms.size();
	ofile << "p cnf " << cnf.var_num << " " << clauses_num << endl;
	for (auto &clause: cnf.clauses) {
		ofile << clause << endl;
	}
	for (auto &clause: clauses_cms) {
		ofile << clause << endl;
	}
	for (auto &clause: clauses_partial_ls) {
		ofile << clause << endl;
	}
	ofile.close();

	// Construct an output file name and create it by opening:
	string local_out_file_name = "./out_cms" + to_string(wu.cms_index) + 
		"_" + "partial-ls" + to_string(wu.partial_ls_index);
	fstream local_out_file;
	local_out_file.open(local_out_file_name, ios_base::out);

	// Run an AllSAT solver once:
	enumerate_by_allsat_solver(solver_name, cur_cnf_name, local_out_file_name);

	vector<vector<int>> sat_assignments = 
		parse_sat_assignments(local_out_file_name, true, ls_order);

	/*
	set<string> first_ls_set;
	for (auto &s_a : sat_assignments) {
		// Typicaly there are N^3 variables, but clasp can produce less:
		if (s_a.size() != 3*pow(ls_order,3)) continue;
		first_ls_set.insert(first_ls_from_sat(s_a, ls_order));
	}
	*/

	// Remove the formed CNF:
	string system_str = "rm " + cur_cnf_name;
	string res_str = exec(system_str);

	wu_result res = UNSAT;
	if (sat_assignments.size() > 0) {
		res = SAT;
	} else {
		// Remove the output file if no solutions:
		system_str = "rm " + local_out_file_name;
		res_str = exec(system_str);
	}

	return res;
}

// Execute a command in a UNIX-based OS:
string exec(const string cmd_str)
{
	assert(cmd_str.size() > 1);
	string result = "";
	char* cmd = new char[cmd_str.size() + 1];
	for (unsigned i = 0; i < cmd_str.size(); i++)
		cmd[i] = cmd_str[i];
	cmd[cmd_str.size()] = '\0';
	FILE* pipe = popen(cmd, "r");
	delete[] cmd;
	if (!pipe) return "ERROR";
	char buffer[128];
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
	pclose(pipe);
	return result;
}


// Parse satysfying assignments from a solver's output:
vector<vector<int>> parse_sat_assignments(
	const string output_str,
	const bool is_file_output,
	const unsigned ls_order
	)
{
	vector<string> lines;
	string line;

	// If an output file's name is given:
	if (is_file_output) {
		ifstream output_file(output_str, ios_base::in);
		while (getline(output_file, line)) {
			lines.push_back(line);
		}
		output_file.close();
	}
	// If an output string is given:
	else {
		if (output_str.find("\n") == string::npos) {
			lines.push_back(output_str);
		}
		else {
			stringstream sstream(output_str);
			while(getline(sstream, line,'\n')) lines.push_back(line);
		}
	}

	assert(lines.size() > 0);

	bool is_started_sol = false;
	vector<vector<int>> all_assignments;
	vector<int> cur_assignment;

	// Process all lines:
	for (auto &line : lines) {
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
					all_assignments.push_back(cur_assignment);
				}
		}
	}

	return all_assignments;
}

// Given a SAT assignment, form the first LS from the found MOLS:
string first_ls_from_sat(const vector<int> sat_assignment,
												 const unsigned ls_order) 
{
	assert(sat_assignment.size() > 0);
	assert(ls_order > 0);
	assert(sat_assignment.size() == 3*pow(ls_order, 3));
	string ls_str = "";
	vector<int> cell_literals;
	for (auto &lit : sat_assignment) {
		// Literals of a current cell are collected:
		cell_literals.push_back(lit);
		int val = -1;
		unsigned neg_lit_num = 0;
		if (cell_literals.size() == ls_order) {
			for (unsigned j=0; j < ls_order; ++j) {
				if (cell_literals[j] > 0) {
					assert(val == -1);
					val = (int)j;
				}
				else {
					neg_lit_num += 1;
				}
			}
			assert(neg_lit_num == ls_order - 1);
			assert(val >= 0 and (unsigned)val <= ls_order - 1);
			ls_str += to_string(val);
			// Read only the first LS:
			if (ls_str.size() == pow(ls_order, 2)) break;
			cell_literals.clear();
		}
	}
	assert(ls_str.size() == ls_order * ls_order);
	return ls_str;
}


// Enumerate all solutions by an AllSAT solver and write them to a file:
void enumerate_by_allsat_solver(string solver_name,
															  const string cur_cnf_name,
															  const string local_out_file_name)
{
	// Find all satisfying assignments of the formed CNF via a SAT solver:
	string solver_params = "";

	// Parse clasp's parameters:
	if (solver_name.find("clasp") != string::npos) {
		string clasp_config_str = "auto";
		string clasp_enum_str = "auto";
		size_t first = solver_name.find("-");
		size_t last = solver_name.find_last_of("-");
		if (first != string::npos and last != string::npos) {
			clasp_config_str = solver_name.substr(first+1, last-first-1);
			clasp_enum_str = solver_name.substr(last+1, solver_name.size()-last-1);
			// Cut the solver name to get an executable clasp name:
			solver_name = solver_name.substr(0, first);
		}
		solver_params = "--configuration=" + clasp_config_str +
				" --enum-mode=" + clasp_enum_str;
	}

	// Enable the solver's enumeration mode:
	if (solver_name.find("cryptominisat") != string::npos) {
		solver_params = "--maxsol 10000000";
	}
	else if (solver_name.find("picosat") != string::npos) {
		solver_params = "--all";
	}
	// Clasp with default enumeration settings:
	else if (solver_name.find("clasp") != string::npos) {
			solver_params += " --models=0";
	}

	// Run the solver:
	string system_str = solver_name + " " + solver_params + " " +
	cur_cnf_name + " > " + local_out_file_name;
	//cout << system_str << endl;
	string res_str = exec(system_str);
}
