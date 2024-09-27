// Created on: 22 Feb 2024
// Author: Oleg Zaikin
// E-mail: zaikin.icc@gmail.com
//
// For a given CNF (that encodes the search for two orthogonal Latin squares),
// a set of cells mapping schemas (CMS), and a set of partial Latin squares, 
// for the Cartesian product of CMSs and partial Latin squares, add to the CNF
// the CMS constraints, substitute the partial Latin square, and enumerate
// all solutions via a SAT solver.
//
// Example:
//   ./solve_sat_cms cryptominisat5 p.cnf cms.txt partls.txt --allsat -cpunum=2
//=============================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <limits>
#include <thread>
#include <set>
#include <math.h>

#include <omp.h>

std::string program = "solve_sat_cms";
std::string version = "0.1.2";

#define time_point_t std::chrono::time_point<std::chrono::system_clock>
#define cms_t std::vector<std::vector<unsigned>>
#define partial_ls_t std::vector<std::vector<std::string>>

std::vector<cms_t> cms_arr;
std::vector<partial_ls_t> partial_ls_arr;

enum status{ NOT_STARTED = -1, IN_PROGRESS = 0, PROCESSED = 1};
enum result{ UNKNOWN = 0, UNSAT = 1, SAT = 2, INTERR = 3 };

struct workunit {
	int id;
	int cms_index;
	int partial_ls_index;
	status stts;
	result rslt;
	double runtime;
	workunit() : id(-1), cms_index(-1), partial_ls_index(-1), stts(NOT_STARTED),
	  rslt(UNKNOWN), runtime(-1) {};
};

struct CNF {
	long long int var_num;
	long long int clause_num;
	std::vector<std::string> clauses;
	CNF() : var_num(0), clause_num(0), clauses() {}
	CNF(std::string cnf_name) : var_num(0), clause_num(0), clauses() {
		read(cnf_name);
	}
	void read(std::string cnf_name) {
		std::ifstream cnf_file(cnf_name, std::ios_base::in);
		std::string str;
		while (getline(cnf_file, str)) {
			if (str.size() == 0 or str[0] == 'p' or str[0] == 'c')
				continue;
			clauses.push_back(str);
			clause_num++;
			std::stringstream sstream;
			sstream << str;
			long long int ival;
			while (sstream >> ival)	var_num = std::max(llabs(ival), var_num);
		}
		cnf_file.close();
	}
	void print() {
		std::cout << "var_num : " << var_num << std::endl;
		std::cout << "clause_num : " << clause_num << std::endl;
	};
};

std::string exec(const std::string cmd_str);
std::vector<cms_t> read_cms(const std::string cms_file_name);
std::vector<partial_ls_t> read_partial_ls(
	const std::string partial_ls_file_name,
  const unsigned ls_order);
unsigned cnf_var_num(const unsigned ls_order, const unsigned ls_index,
                     const unsigned row_index, const unsigned col_index,
										 const unsigned cell_val);
std::string str_after_prefix(const std::string str, const std::string prefix);
std::vector<std::string> find_all_mols(std::string solver_name,
																			 const bool isAllSATsolver,
																			 const CNF cnf,
																			 workunit &wu,
																			 const time_point_t program_start);
std::vector<std::vector<int>> parse_sat_assignments(
	const std::string output_str, const bool is_file_output,
	const unsigned ls_order);
std::string first_ls_from_sat(const std::vector<int> sat_assignment,
															const unsigned ls_order);
void enumerate_by_sat_solver(const std::string solver_name,
														 const std::string cur_cnf_name,
														 const unsigned ls_order,
														 const std::vector<int> partial_ls_literals,
														 const std::string local_out_file_name);
void enumerate_by_allsat_solver(std::string solver_name,
														    const std::string cur_cnf_name,
														    const std::string local_out_file_name);							 

void print_usage() {
	std::cout << "Usage : ./" << program
						<< " solver CNF CMS partial-LS [Options]\n";
	std::cout << " solver     	  : SAT solver\n";
  std::cout << " CNF        	  : file with CNF that encodes 2MOLS\n";
  std::cout << " CMS        	  : file with list of CMSs\n";
  std::cout << " partial-LS 	  : file with list of partial LSs\n";
	std::cout << "    Options\n";
	std::cout << "    -cpunum=<int> : (default = all cores) CPU cores\n";
	std::cout << "    --allsat      : whether an AllSAT solver is given\n"; 
}

void print_version() {
	std::cout << program << " of version " << version << std::endl;
}

int main(const int argc, const char *argv[]) {
  std::vector<std::string> str_argv;
	for (int i=0; i < argc; ++i) str_argv.push_back(argv[i]);
	assert(str_argv.size() == argc);

	if (argc == 2 and str_argv[1] == "-v") {
		print_version();
		std::exit(EXIT_SUCCESS);
	}

	if (argc == 2 and str_argv[1] == "-h") {
		print_usage();
		std::exit(EXIT_SUCCESS);
	}

  if (argc < 5) {
		print_usage();
		std::exit(EXIT_FAILURE);
	}

	std::string solver_name	         = str_argv[1];
	std::string cnf_file_name        = str_argv[2];
	std::string cms_file_name        = str_argv[3];
	std::string partial_ls_file_name = str_argv[4];
	std::cout << "solver_name          : " << solver_name   << std::endl;
	std::cout << "cnf_file_name        : " << cnf_file_name << std::endl;
	std::cout << "cms_file_name        : " << cms_file_name << std::endl;
	std::cout << "partial_ls_file_name : " << partial_ls_file_name << std::endl;

	unsigned cpu_num = 0;
	bool isAllSATsolver = false;
	if (str_argv.size() > 5) {
		for (unsigned i=5; i < str_argv.size(); ++i) {
				std::string s = str_after_prefix(str_argv[i], "-cpunum=");
				if (s != "") std::istringstream(s) >> cpu_num;
				if (str_argv[i] == "--allsat") isAllSATsolver = true;
		}
	}
	std::cout << "cpu_num              : " << cpu_num << std::endl;
	std::cout << "isAllSATsolver       : " << isAllSATsolver << std::endl;

	// If no CPU num is given, use all threads:
	if (!cpu_num) {
		cpu_num = std::thread::hardware_concurrency();
		std::cout << cpu_num << " CPU cores are detected" << std::endl;
	}

	omp_set_num_threads(cpu_num);
	std::cout << cpu_num << " threads are used" << std::endl;

	// Read a CNF and print its statistics:
	CNF cnf(cnf_file_name);
	cnf.print();

	std::cout << std::endl;

	// Read all CMSs from file:
	cms_arr = read_cms(cms_file_name);
	assert(cms_arr.size() > 0);
	std::cout << cms_arr.size() << " CMSs were read" << std::endl;
	std::cout << "The first CMS is :" << std::endl;
	for (auto x : cms_arr[0]) {
		for (auto y : x) {
			std::cout << y << " ";
		}
		std::cout << std::endl;
	}
	unsigned ls_order = cms_arr[0].size();
	std::cout << std::endl;
	std::cout << "LS order : " << ls_order << " was detected in CMS" << std::endl;

	std::cout << std::endl;
	partial_ls_arr = read_partial_ls(partial_ls_file_name, ls_order);
	std::cout << partial_ls_arr.size() << 
		" partial Latin squares were read" << std::endl;
	std::cout << "The first one is" << std::endl;
	for (auto x : partial_ls_arr[0]) {
		for (auto y : x) {
			std::cout << y << " ";
		}
		std::cout << std::endl;
	}

	std::vector<workunit> wu_arr;
	unsigned wu_id = 0;
	for (unsigned i=0; i<cms_arr.size(); ++i) {
		for (unsigned j=0; j<partial_ls_arr.size(); ++j) {
			workunit wu;
			wu.id = wu_id;
			wu.cms_index = i;
			wu.partial_ls_index = j;
			wu.stts = NOT_STARTED;
			wu.rslt = UNKNOWN;
			wu.runtime = -1;
			wu_arr.push_back(wu);
			//
			wu_id++;
		}
	}
	assert(wu_arr.size() > 0);
	assert(wu_arr.size() == cms_arr.size() * partial_ls_arr.size());
	std::cout << wu_arr.size() << " workunits were generated" << std::endl;

	const time_point_t program_start = std::chrono::system_clock::now();
	unsigned solved_wus = 0;
	std::set<std::string> all_first_lss;

	// Process all workunits in parallel:
	#pragma omp parallel for schedule(dynamic, 1)
	for (auto &wu : wu_arr) {
		std::vector<std::string> wu_first_lss = 
			find_all_mols(solver_name, isAllSATsolver, cnf, wu, program_start);
		for (auto &ls : wu_first_lss) all_first_lss.insert(ls);
		solved_wus++;
		if (ls_order > 8 or solved_wus % 1000 == 0) {
			std::cout << solved_wus << " workunits out of " << wu_arr.size()
								<< " are solved\n";
		}
	}
	std::cout << all_first_lss.size() << " MOLS were found\n";

	const time_point_t program_end = std::chrono::system_clock::now();
	std::cout << "Elapsed : "
	<< std::chrono::duration_cast<std::chrono::seconds>(program_end - program_start).count()
	<< " seconds" << std::endl;

	size_t pos = solver_name.find("/");
	std::string base_solver_name = solver_name.substr(pos+1, solver_name.size() - pos);
	std::string log_file_name = "log_" + base_solver_name + "_n" + std::to_string(ls_order);
	std::cout << "Writing log to file " << log_file_name << std::endl;
	std::ofstream log_file(log_file_name, std::ios_base::out);
	log_file << "Input parameters:" << std::endl;
	for (auto &str : str_argv) log_file << str << std::endl;
	log_file << std::endl;
	log_file << solved_wus << " workunits out of " << wu_arr.size()
								<< " are solved" << std::endl;
	log_file << all_first_lss.size() << " MOLS were found" << std::endl;
	log_file << "Elapsed : "
	<< std::chrono::duration_cast<std::chrono::seconds>(program_end - program_start).count()
	<< " seconds" << std::endl;
	log_file.close();

	std::string esodls_file_name = "esodls_n" + std::to_string(ls_order) + ".txt";
	std::ofstream esodls_file(esodls_file_name, std::ios_base::out);
	std::cout << "Writing found ESODLS to file " << esodls_file_name << std::endl;
	for (auto &ls : all_first_lss) {
		esodls_file << ls << "\r\n";
	}
	esodls_file.close();

	return 0;
}

// Find substring after a given prefix:
std::string str_after_prefix(const std::string str, const std::string prefix) {
    size_t pos = str.find( prefix );
    if ( pos != std::string::npos )
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

// Enumerate all solutions by an AllSAT solver and write them to a file:
void enumerate_by_allsat_solver(std::string solver_name,
															  const std::string cur_cnf_name,
															  const std::string local_out_file_name)
{
	// Find all satisfying assignments of the formed CNF via a SAT solver:
	std::string solver_params = "";

	// Parse clasp's parameters:
	if (solver_name.find("clasp") != std::string::npos) {
		std::string clasp_config_str = "auto";
		std::string clasp_enum_str = "auto";
		std::size_t first = solver_name.find("-");
		std::size_t last = solver_name.find_last_of("-");
		if (first != std::string::npos and last != std::string::npos) {
			clasp_config_str = solver_name.substr(first+1, last-first-1);
			clasp_enum_str = solver_name.substr(last+1, solver_name.size()-last-1);
			// Cut the solver name to get an executable clasp name:
			solver_name = solver_name.substr(0, first);
		}
		solver_params = "--configuration=" + clasp_config_str +
				" --enum-mode=" + clasp_enum_str;
	}

	// Enable the solver's enumeration mode:
	if (solver_name.find("cryptominisat") != std::string::npos) {
		solver_params = "--maxsol 10000000";
	}
	else if (solver_name.find("picosat") != std::string::npos) {
		solver_params = "--all";
	}
	// Clasp with default enumeration settings:
	else if (solver_name.find("clasp") != std::string::npos) {
			solver_params += " --models=0";
	}

	// Run the solver:
	std::string system_str = solver_name + " " + solver_params + " " +
	cur_cnf_name + " > " + local_out_file_name;
	//std::cout << system_str << std::endl;
	std::string res_str = exec(system_str);
}

// Enumerate all solutions by a SAT solver and write them to a file:
void enumerate_by_sat_solver(const std::string solver_name,
														 const std::string cur_cnf_name,
														 const unsigned ls_order,
														 const std::vector<int> partial_ls_literals,
														 const std::string local_out_file_name)
{
	CNF cnf(cur_cnf_name);
	std::vector<std::vector<int>> all_sat_assignments; 
	bool UNSAT = false;
	std::vector<std::string> blocking_clauses;
	//std::cout << solver_name << " on CNF " << cur_cnf_name << std::endl;
	// Run the SAT solver, add the blocking clause, and run again:
	while (not UNSAT) {
		std::string system_str = solver_name + " " + cur_cnf_name;
		//std::cout << system_str << std::endl;
		std::string solver_output_str = exec(system_str);
		assert(solver_output_str.size() > 0);
		std::vector<std::vector<int>> sat_assignments = 
			parse_sat_assignments(solver_output_str, false, ls_order);
		assert(sat_assignments.size() <= 1);
		if (sat_assignments.size() == 0 ) { 
			UNSAT = true;
		}
		else {
			all_sat_assignments.push_back(sat_assignments[0]);
			// Form the blocking clause:
			std::string blocking_clause;
			for (auto &lit : sat_assignments[0]) {
				//blocking_clause += std::to_string(-lit) + " "; // all literals
				//if (lit > 0 and lit <= (int)pow(ls_order, 3))  // positive first-ls literals
				// positive first-ls literals which are not in partial ls:
				if (lit > 0 and lit <= (int)pow(ls_order, 3) and 
					std::find(partial_ls_literals.begin(), partial_ls_literals.end(), lit) == partial_ls_literals.end())
					blocking_clause += std::to_string(-lit) + " "; 
			}
			blocking_clause += "0";
			//std::cout << blocking_clause << std::endl;
			assert(blocking_clause.size() > 0);
			blocking_clauses.push_back(blocking_clause);
			// Write the CNF from the scratch:
			std::ofstream cur_cnf_file(cur_cnf_name, std::ios_base::out);
			cur_cnf_file << "p cnf " << cnf.var_num << " "
				           << cnf.clause_num + blocking_clauses.size() << std::endl;
			for (auto &c : cnf.clauses) cur_cnf_file << c << std::endl;
			// Add all blocking clauses formed so far:
			for (auto &c : blocking_clauses) cur_cnf_file << c << std::endl;
			cur_cnf_file.close();
		}
	}

	// Write all solutions to a file:
	std::ofstream local_out_file(local_out_file_name, std::ios_base::out);
	if (all_sat_assignments.size() > 0) {
		for (auto &s_a : all_sat_assignments) {
			local_out_file << "s SATISFIABLE\n";
			local_out_file << "v ";
			for (auto &lit : s_a) {
				local_out_file << lit << " ";
			}
			local_out_file << "0" << std::endl;
		}
	}
	else {
		local_out_file << "s UNSATISFIABLE\n";
	}
	local_out_file.close();
}

// Find all MOLS in a CNF with added CMS- and partial-LS-constraints:
std::vector<std::string> find_all_mols(std::string solver_name,
																			 const bool isAllSATsolver,
																			 const CNF cnf,
																			 workunit &wu,
																			 const time_point_t program_start)
{
	unsigned ls_order = partial_ls_arr[wu.partial_ls_index].size();

	// CMS constraints:
	std::vector<std::string> clauses_cms;
  for (unsigned i=0; i < ls_order; ++i) {
		for (unsigned j=0; j < ls_order; ++j) {
			std::vector<unsigned> first_ls_cell_vars;
			// CNF variables for values in cell (i,j) of LS 1
			for (unsigned k=0; k < ls_order; ++k) {
				first_ls_cell_vars.push_back(cnf_var_num(ls_order, 0, i, j, k));
			}
			assert(first_ls_cell_vars.size() == ls_order);
			// cell (i, j) in LS 1 == cell (i2, j2) in LS 2:
			unsigned i2 = cms_arr[wu.cms_index][i][j] / ls_order;
      unsigned j2 = cms_arr[wu.cms_index][i][j] % ls_order;
			std::vector<unsigned> second_ls_cell_vars;
			// CNF variables for values in cell (i2,j2) of LS 2
			for (unsigned k=0; k < ls_order; ++k) {
				second_ls_cell_vars.push_back(cnf_var_num(ls_order, 1, i2, j2, k));
			}
			assert(second_ls_cell_vars.size() == ls_order);
			// Equality-clauses for each pair of values:
			for (unsigned k=0; k < ls_order; ++k) {
				std::string str_val1 = std::to_string(first_ls_cell_vars[k]);
				std::string str_val2 = std::to_string(second_ls_cell_vars[k]);
				clauses_cms.push_back(      str_val1 + " -" + str_val2 + " 0");
				clauses_cms.push_back("-" + str_val1 + " "  + str_val2 + " 0");
			}
		}
	}
	assert(clauses_cms.size() == pow(ls_order,3)*2);

	// Partial-LS constraints for the first LS:
	std::vector<int> partial_ls_literals;
	std::vector<std::string> clauses_partial_ls;
	for (unsigned i=0; i < ls_order; ++i) {
		for (unsigned j=0; j < ls_order; ++j) {
			std::string val_str = partial_ls_arr[wu.partial_ls_index][i][j];
			if (val_str == "-") continue;
			unsigned val_uint;
			std::istringstream(val_str) >> val_uint;
			unsigned var = cnf_var_num(ls_order, 0, i, j, val_uint);
			clauses_partial_ls.push_back(std::to_string(var) + " 0");
			partial_ls_literals.push_back(var);
		}
	}

	std::string cur_cnf_name = "./cms" + std::to_string(wu.cms_index) + "_" +
		"partial-ls" + std::to_string(wu.partial_ls_index) + ".cnf";

	// Write all needed clauses to a file:
	std::ofstream ofile(cur_cnf_name, std::ios_base::out);
	unsigned clauses_num = cnf.clauses.size() + clauses_partial_ls.size() +
		clauses_cms.size();
	ofile << "p cnf " << cnf.var_num << " " << clauses_num << std::endl;
	for (auto &clause: cnf.clauses) {
		ofile << clause << std::endl;
	}
	for (auto &clause: clauses_cms) {
		ofile << clause << std::endl;
	}
	for (auto &clause: clauses_partial_ls) {
		ofile << clause << std::endl;
	}
	ofile.close();

	// Construct an output file name and create it by opening:
	std::string local_out_file_name = "./out_cms" +
		std::to_string(wu.cms_index) + "_" + "partial-ls" + 
		std::to_string(wu.partial_ls_index);
	std::fstream local_out_file;
	local_out_file.open(local_out_file_name, std::ios_base::out);

	const time_point_t t1 = std::chrono::system_clock::now();
	
	// Run an AllSAT solver once or a SAT solver iteratively:
	if (isAllSATsolver) {
		enumerate_by_allsat_solver(solver_name, cur_cnf_name, local_out_file_name);
	}
	else {
		enumerate_by_sat_solver(solver_name, cur_cnf_name, ls_order, 
			partial_ls_literals, local_out_file_name);
	}

	std::vector<std::vector<int>> sat_assignments = 
		parse_sat_assignments(local_out_file_name, true, ls_order);

	if (sat_assignments.size() > 0) {
		const time_point_t t2 = std::chrono::system_clock::now();
		//std::cout << solver_name << " on CNF " << cur_cnf_name <<
		std::cout << cur_cnf_name << " SAT : "
		<< std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
		<< " seconds" << std::endl;
	}

	// Remove the formed CNF:
	std::string system_str = "rm " + cur_cnf_name;
	std::string res_str = exec(system_str);

	// Remove the output file:
	system_str = "rm " + local_out_file_name;
	res_str = exec(system_str);

	std::vector<std::string> first_ls_arr;
	for (auto &s_a : sat_assignments) {
		// Typicaly there are N^3 variables, but clasp can produce less:
		if (s_a.size() != 3*pow(ls_order,3)) continue;
		first_ls_arr.push_back(first_ls_from_sat(s_a, ls_order));
	}

	return first_ls_arr;
}

// Execute a system command in a Linux OS:
std::string exec(const std::string cmd_str) {
	char* cmd = new char[cmd_str.size() + 1];
	for (unsigned i = 0; i < cmd_str.size(); i++)
		cmd[i] = cmd_str[i];
	cmd[cmd_str.size()] = '\0';
	FILE* pipe = popen(cmd, "r");
	delete[] cmd;
	if (!pipe) return "ERROR";
	char buffer[128];
	std::string result = "";
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
	pclose(pipe);
	return result;
}

// Parse satysfying assignments from a solver's output:
std::vector<std::vector<int>> parse_sat_assignments(
	const std::string output_str,
	const bool is_file_output,
	const unsigned ls_order
	)
{
	std::vector<std::string> lines;
	std::string line;

	// If an output file's name is given:
	if (is_file_output) {
		std::ifstream output_file(output_str, std::ios_base::in);
		while (getline(output_file, line)) {
			lines.push_back(line);
		}
		output_file.close();
	}
	// If an output string is given:
	else {
		if (output_str.find("\n") == std::string::npos) {
			lines.push_back(output_str);
		}
		else {
			std::stringstream sstream(output_str);
			while(std::getline(sstream, line,'\n')) lines.push_back(line);
		}
	}

	assert(lines.size() > 0);

	bool is_started_sol = false;
	std::vector<std::vector<int>> all_assignments;
	std::vector<int> cur_assignment;

	// Process all lines:
	for (auto &line : lines) {
		if (line.size() < 2) continue;
		if ( (line.find("s SATISFIABLE") != std::string::npos) or
		     (line.find("c Answer: ") != std::string::npos) )
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
				std::stringstream sstream(line);
				std::vector<int> literals;
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
		// cadical-allsat's format:
		if (line.find("c New solution:") != std::string::npos) {
			// Delete the first 3 words of total length 16:
			line.erase(0, 16);
			// Read literals
			std::stringstream sstream(line);
			std::vector<int> literals;
			int ival;
			assert(cur_assignment.empty());
			while (sstream >> ival) cur_assignment.push_back(ival);
			assert(cur_assignment[cur_assignment.size() - 1] == 0);
			cur_assignment.pop_back();
			all_assignments.push_back(cur_assignment);
			cur_assignment.clear();
		}
	}

	return all_assignments;
}

// Given a SAT assignment, form the first LS from the found MOLS:
std::string first_ls_from_sat(const std::vector<int> sat_assignment,
															const unsigned ls_order) 
{
	assert(sat_assignment.size() > 0);
	assert(ls_order > 0);
	assert(sat_assignment.size() == 3*pow(ls_order, 3));
	std::string ls_str = "";
	std::vector<int> cell_literals;
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
			ls_str += std::to_string(val);
			// Read only the first LS:
			if (ls_str.size() == pow(ls_order, 2)) break;
			cell_literals.clear();
		}
	}
	assert(ls_str.size() == ls_order * ls_order);
	return ls_str;
}

std::vector<cms_t> read_cms(const std::string cms_file_name) {
	std::vector<cms_t> cms_arr;
	std::ifstream cms_file(cms_file_name, std::ios_base::in);
	assert(cms_file.is_open());
	std::string str;
	unsigned ls_order = 0;
	cms_t cms;
	while (getline(cms_file, str)) {
		if ( (str.size() < 2) or (str[0] == '#') ) continue;
		if ( (str.find("Loading") != std::string::npos) or
			   (str.find("ESODLS") != std::string::npos)
				) continue;
		std::vector<unsigned> row;
		std::stringstream sstream;
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

std::vector<partial_ls_t> read_partial_ls( const std::string partial_ls_file_name,
  const unsigned ls_order)
{
	std::vector<partial_ls_t> partial_ls_arr;
	assert(ls_order);

	std::ifstream partial_ls_file(partial_ls_file_name, std::ios_base::in);
	assert(partial_ls_file.is_open());
	std::string str;
	partial_ls_t partial_ls;
	while (getline(partial_ls_file, str)) {
		if ( (str.size() < 2) or 
		     (str.find("filling") != std::string::npos) 
			 )
				continue;
		std::vector<std::string> row;
		std::stringstream sstream;
		sstream << str;
		std::string word;
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
