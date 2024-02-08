# Created on: 8 Feb 2024
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# Given a path with files where a SAT solver's results are stored, collect all
# satisfying assignments and build pairs of MOLS based on them. 
#=============================================================================

import sys
import os

script = "collect_mols_from_sat.py"
version = "0.0.1"

def ls_from_sat(ls_sol : list, ls_order : int):
	assert(len(ls_sol) > 0)
	assert(ls_order > 0)
	ls = []
	cell_literals = []
	for lit in ls_sol:
		cell_literals.append(lit)
		# Literals of a current cell are collected:
		val = -1
		neg_lit_num = 0
		if len(cell_literals) == ls_order:
			for j in range(ls_order):
				if int(cell_literals[j]) > 0:
					val = j
				else:
					neg_lit_num += 1
			assert(neg_lit_num == ls_order - 1)
			assert(val >= 0 and val <= ls_order - 1)
			ls.append(val)
			cell_literals = []
	assert(len(ls) == ls_order * ls_order)
	return ls

def parse_sat_assignments(fname : str):
	is_started_sol = False
	sat_assignments = []
	cur_assignment = []
	with open(fname, 'r') as ifile:
		lines = ifile.read().splitlines()
		for line in lines:
			if len(line) < 2:
				continue
			if 's SATISFIABLE' in line:
				is_started_sol = True
				cur_assignment = []
				continue
			if is_started_sol:
				print(fname)
				print(line)
				assert(line[:2] == 'v ')
				# If the last line of a solution is found, stop reading:
				if line[-2:] == ' 0':
					is_started_sol = False
					sat_assignments.append(cur_assignment)
					# Cut '0' at the end:
					line = line[:-1]
				# Cut 'v' at the beginning:
				line = line[1:]
				# Read literals
				literals = line.split()
				assert(len(literals) > 0)
				cur_assignment += literals
	return sat_assignments

### Main function:
print(script + ' of version ' + version + ' is running')
if len(sys.argv) == 2 and sys.argv[1] == '-v':
	exit(1)

if len(sys.argv) < 3:
	print('Usage : ' + script + ' ls-order path-sat-solutions')
	print('  ls-order           : order of LS in MOLS')
	print('  path-sat-solutions : path where SAT solutions are stored')
	exit(1)

ls_order = int(sys.argv[1])
path = sys.argv[2]
print('ls_order : ' + str(ls_order))
print('path     : ' + path)
cur_dir = os.getcwd()
os.chdir(path)
res_filenames = []
for f in os.listdir(path):
    if f.endswith(".out"):
        res_filenames.append(f)
print(str(len(res_filenames)) + ' files with SAT solver results are found')

one_ls_var_num = pow(ls_order, 3)
mols_lst = []
for fname in res_filenames:
	sat_assignments = parse_sat_assignments(fname)
	if len(sat_assignments) == 0:
		continue
	for assign in sat_assignments:
		ls1_sol = assign[:one_ls_var_num]
		ls2_sol = assign[one_ls_var_num:2*one_ls_var_num]
		ls1 = ls_from_sat(ls1_sol, ls_order)
		ls2 = ls_from_sat(ls2_sol, ls_order)
		mols_lst.append([ls1, ls2])

print(str(len(mols_lst)) + ' mols in total')

# Go back to the current directory:
os.chdir(cur_dir)
esodls_fname = 'esodls_n' + str(ls_order)
print('Writing to file ' + esodls_fname)
with open(esodls_fname, 'w') as ofile:
	for mols in mols_lst:
		for i in range(ls_order):
			for j in range(ls_order):
				s = str(mols[0][i*ls_order + j]) + ' '
				ofile.write(s)
			ofile.write('\n')
		ofile.write('\n')
