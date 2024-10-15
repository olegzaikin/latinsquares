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
version = "0.0.9"

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

def merge_assign_lines(lines : list):
	for i in range(len(lines) - 1):
		if len(lines[i]) < 1 or len(lines[i+1]) < 1:
			continue
		# If i-th string is not the last line with a solution,
		#   and (i+1)-th string does not start with 'v', merge them
		if lines[i][0] == 'v' and (len(lines[i]) == 1 or lines[i][-2:] != ' 0')\
			and lines[i+1][0] != 'v':
			lines[i] += lines[i+1]
			lines[i+1] = ''
	return lines

def ascending_order(literals : list):
	vars = [abs(int(l)) for l in literals]
	return vars == sorted(vars)

def parse_sat_assignments(fname : str, ls_order : int):
	is_started_sol = False
	sat_assignments = []
	cur_assignment = []
	cms_indices = []
	with open(fname, 'r') as ifile:
		lines = ifile.read().splitlines()
		# Process divided lines with sat assignments:
		lines = merge_assign_lines(lines)
		cms_index = -1
		for line in lines:
			if len(line) < 2:
				continue
			#c Reading from 2mols_n7_diag_cms=000_filling=0.cnf
			if 'c Reading from ' in line and '_cms=' in line and '_' in line:
				cms_index = int(line.split('_cms=')[1].split('_')[0])
			if 's SATISFIABLE' in line or 'c Answer: ' in line:
				is_started_sol = True
				cur_assignment = []
				continue
			if is_started_sol:
				# If clasp solver, no solution after after 's SATISFIABLE':
				if line[0] != 'v' or line[1] != ' ':
					continue
				# Cut 'v' at the beginning:
				line = line[1:]
				# Read literals
				literals = line.split()
				assert(ascending_order(literals) or literals[-1] == '0')
				assert(len(literals) > 0)
				cur_assignment += literals
				# If the last line of a solution is found, stop reading:
				if literals[-1] == '0':
					is_started_sol = False
					# Cut 0 at the end:
					assert(cur_assignment[-1] == '0')
					cur_assignment = cur_assignment[:-1]
					# The number of variables must be 3*N^3:
					#print(str(len(cur_assignment)))
					assert(len(cur_assignment) == 3*pow(ls_order, 3))
					sat_assignments.append(cur_assignment)
					if cms_index >= 0:
						cms_indices.append(cms_index)
	return sat_assignments, cms_indices

def ls_as_str(ls : list, ls_order : int):
	res_str = ''
	for i in range(ls_order):
		for j in range(ls_order):
			res_str += str(ls[i*ls_order + j])
			if j < ls_order-1:
				res_str += ' '
		res_str += '\n'
	res_str += '\n'
	return res_str

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
for dirpath,_,filenames in os.walk(path):
	for f in filenames:
		print(f)
		if f.startswith("!sat_out") or f.startswith("out_cms"):
			res_filenames.append(os.path.abspath(os.path.join(dirpath, f)))
#for f in os.listdir(path):
    #if f.startswith("!sat_out"):
        #res_filenames.append(f)
print(str(len(res_filenames)) + ' files with SAT solver results are found')

one_ls_var_num = pow(ls_order, 3)
mols_lst = []
k = 0
cms_indices = []
for fname in res_filenames:
	if k % 100 == 0 and k > 0:
		print(str(k) + ' files out of ' + str(len(res_filenames)) + ' are processed')
	#print(fname)
	sat_assignments, cms_indices = parse_sat_assignments(fname, ls_order)
	k += 1
	if len(sat_assignments) == 0:
		continue
	for assign in sat_assignments:
		ls1_sol = assign[:one_ls_var_num]
		ls2_sol = assign[one_ls_var_num:2*one_ls_var_num]
		ls1 = ls_from_sat(ls1_sol, ls_order)
		ls2 = ls_from_sat(ls2_sol, ls_order)
		mols_lst.append([ls1, ls2])

print(str(len(mols_lst)) + ' mols in total')
print('cms_indices :')
print(cms_indices)

esodls_str = set()
mols_lst_str_lst = []
sodls_lst = []
esodls_not_sodls_lst = []
for i in range(len(mols_lst)):
	s = ''
	for val in mols_lst[i][0]:
		s += str(val)
	esodls_str.add(s)
	#s = ''
	#for val in mols[1]:
	#	s += str(val)
	#esodls_str.add(s)
	# All pairs of MOLS to string:
	if len(cms_indices) > i:
		assert(cms_indices[i] > 0)
		mols_lst_str_lst.append('cms_index ' + str(cms_indices[i]) + '\n')
	for mols_ind in range(2):
		if mols_lst[i][mols_ind] not in mols_lst:
			mols_lst_str_lst.append(ls_as_str(mols_lst[i][mols_ind], ls_order))
	mols_lst_str_lst.append('\n')
	if len(cms_indices) == 0:
		continue
	# If cms1, then sodls:
	if cms_indices[i] == 1:
		if mols_lst[i][0] not in sodls_lst:
			sodls_lst.append(mols_lst[i][0])

if len(sodls_lst) > 0:
	for i in range(len(mols_lst)):
		if cms_indices[i] > 1 and mols_lst[i][0] not in sodls_lst and mols_lst[i][0] not in esodls_not_sodls_lst:
			esodls_not_sodls_lst.append(mols_lst[i][0])

sorted_esodls_str = sorted(esodls_str)

print(str(len(esodls_str)) + ' esodls in total')

# Go back to the current directory:
os.chdir(cur_dir)
esodls_fname = 'esodls_n' + str(ls_order) + '.txt'
print('Writing to file ' + esodls_fname)
with open(esodls_fname, 'w') as ofile:
	for ls_str in sorted_esodls_str:
		ofile.write(ls_str + '\r\n')

if len(mols_lst_str_lst) > 0:
	mols_fname = 'mols_esodls_n' + str(ls_order) + '.txt'
	print('Writing to file ' + mols_fname)
	with open(mols_fname, 'w') as ofile:
		for s in mols_lst_str_lst:
			ofile.write(s)

if len(sodls_lst) > 0:
	mols_fname = 'mols_sodls_n' + str(ls_order) + '.txt'
	print('Writing to file ' + mols_fname)
	with open(mols_fname, 'w') as ofile:
		for ls in sodls_lst:
			ofile.write(ls_as_str(ls, ls_order))

if len(esodls_not_sodls_lst) > 0:
	mols_fname = 'mols_esodls_not_sodls_n' + str(ls_order) + '.txt'
	print('Writing to file ' + mols_fname)
	with open(mols_fname, 'w') as ofile:
		for ls in esodls_not_sodls_lst:
			ofile.write(ls_as_str(ls, ls_order))
