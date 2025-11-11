# Created on: 7 Nov 2026
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# Encode to SAT the search for transversals of a Latin square of order n.
#
# Example:
# 	python3 ./sat_transversals.py 4 4 --diag
#=============================================================================

import sys
import itertools
import math
import os
import itertools

script = "sat_transversals.py"
version = "0.0.2"

# Return clauses that encode the AtMostOne constraint:
def at_most_one_clauses(vars : list):
	assert(len(vars) > 1)
	assert(len(set(vars)) == len(vars))
	res_clauses = []
	for i in range(len(vars)):
		for j in range(i+1, len(vars)):
			res_clauses.append([-vars[i], -vars[j]])
	return res_clauses

# Return a clause that encodes the AtLeastOne constraint:
def at_least_one_clauses(vars : list):
	res_clauses = [vars]
	return res_clauses

# Get clauses encoding the ExactlyOne constraint for a set of variables:
def exactly_one_clauses(vars : list):
	assert(len(vars) > 1)
	assert(len(set(vars)) == len(vars))
	res_clauses = at_most_one_clauses(vars)
	res_clauses += at_least_one_clauses(vars)
	return res_clauses

# Latin square constraints:
def latin_square_clauses(A : list, is_diag : bool):
	ls_order = len(A)
	assert(ls_order > 0)
	res_clauses = []
	# Constraints on rows, columns, and values are obligatory:
	for i in range(ls_order):
		for j in range(ls_order):
			# Each square's cell contains exactly one value 0..n-1:
			res_clauses += exactly_one_clauses([A[i][j][k] for k in range(ls_order)])
			# Each value occurs exactly once in each row:
			res_clauses += exactly_one_clauses([A[i][k][j] for k in range(ls_order)])
			# Each value occurs exactly once in each column:
			res_clauses += exactly_one_clauses([A[k][i][j] for k in range(ls_order)])
	# Constraints on main diagonal and antidiagonal are optional:
	if is_diag:
		# Main diagonal:
		for i in range(ls_order):
			res_clauses += exactly_one_clauses([A[k][k][i] for k in range(ls_order)])
		# Main antidiagonal:
		for i in range(ls_order):
			res_clauses += exactly_one_clauses([A[ls_order-k-1][k][i] for k in range(ls_order)])
	return res_clauses

def first_row_asending_clauses(A : list):
	ls_order = len(A)
	assert(ls_order > 0)
	res_clauses = []
	for j in range(ls_order):
		res_clauses.append([A[0][j][j]])
	return res_clauses

def known_main_diagonal_clauses(A : list):
	ls_order = len(A)
	assert(ls_order > 0)
	res_clauses = []
	# Main diagonal is in the ascending order 0, 1, ... , n-1
	for i in range(ls_order):
		res_clauses.append([A[i][i][i]])
	perms = list(itertools.permutations([i for i in range(ls_order)]))
	#print(perms)
	possible_main_antidiags = []
	for p in perms:
		is_suit = True
		for j in range(len(p)):
			if p[j] == j or p[j] == ls_order-1-j:
				is_suit = False
				break
		if is_suit:
			possible_main_antidiags.append(p)
	#print(possible_main_antidiags)
	assert(len(possible_main_antidiags) > 0)
	print(str(len(possible_main_antidiags)) + ' possible main antidiagonals')
	# Main antidiagonal is also known.
	# This is the first possible value if the main diagonal is in the ascending order.
	for i in range(ls_order):
		res_clauses.append([A[ls_order-1-i][i][possible_main_antidiags[0][i]]])
	return res_clauses

def transversals_clauses(ls_order : int, ls_variables : list, cpindex_variables : list, transversals_variables : list, is_diag : bool):
	assert(ls_order > 0)
	transversals_num = len(transversals_variables)
	assert(transversals_num > 0)
	assert(len(cpindex_variables) == transversals_num)
	res_clauses = []
	# Each transversal is encoded by N distinct columns' indices (rows indices are in the ascending order):
	A = ls_variables
	for trv_index in range(transversals_num):
		T = transversals_variables[trv_index]
		Z = cpindex_variables[trv_index]
		for j in range(ls_order):
			# Each transversal's cell contains exactly one value (column index) 0..n-1:
			res_clauses += exactly_one_clauses([T[j][k] for k in range(ls_order)])
			# Each column index occurs once among all transversal's cell
			res_clauses += exactly_one_clauses([T[k][j] for k in range(ls_order)])
		for j in range(ls_order):
			# Each cpindex's cell contains exactly one value (column index) 0..n-1:
			res_clauses += exactly_one_clauses([Z[j][k] for k in range(ls_order)])
			# Each column index occurs once among all cpindex's cell
			res_clauses += exactly_one_clauses([Z[k][j] for k in range(ls_order)])
		# Now the transversals must to be connected with the Latin square.
		# Element indexing constraint Z[T_i] = Aij.
		for i in range(ls_order):
			for j in range(ls_order):
				for k in range(ls_order):
					res_clauses.append([-A[i][j][k], -T[i][j], Z[j][k]])
					res_clauses.append([-T[i][j], -Z[j][k], A[i][j][k]])
					res_clauses.append([-A[i][j][k], -Z[j][k], T[i][j]])
		if is_diag:
			# Exactly one entry from the main diagonal:
			res_clauses += exactly_one_clauses([T[j][j] for j in range(ls_order)])
			# Exactly one entry from the main antidiagonal:
			res_clauses += exactly_one_clauses([T[j][ls_order - 1 - j] for j in range(ls_order)])
	return res_clauses

def different_transversals_clauses(ls_order : int, T1 : list, T2 : list, xor_variables : list):
	trv_len = len(T1)
	assert(trv_len > 0)
	assert(trv_len == ls_order)
	assert(len(T2) == trv_len and len(xor_variables) == trv_len)
	res_clauses = []
	xor_clause = []
	for i in range(ls_order):
		for j in range(ls_order):
			res_clauses.append([-T1[i][j], -T2[i][j], -xor_variables[i][j]])
			res_clauses.append([-T1[i][j], T2[i][j], xor_variables[i][j]])
			res_clauses.append([T1[i][j], -T2[i][j], xor_variables[i][j]])
			res_clauses.append([T1[i][j], T2[i][j], -xor_variables[i][j]])
			xor_clause.append(xor_variables[i][j])
	res_clauses.append(xor_clause)
	return res_clauses

def ascending_order(literals : list):
	vars = [abs(int(l)) for l in literals]
	return vars == sorted(vars)

def parse_sat_assignments(fname : str, ls_order : int):
	is_started_sol = False
	sat_assignment = []
	with open(fname, 'r') as ifile:
		lines = ifile.read().splitlines()
		# Process divided lines with sat assignments:
		for line in lines:
			if len(line) < 2:
				continue
			if 's UNSATISFIABLE' in line:
				break
			if 's SATISFIABLE' in line or 'c Answer: ' in line:
				is_started_sol = True
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
				sat_assignment += literals
				# If the last line of a solution is found, stop reading:
				if literals[-1] == '0':
					is_started_sol = False
					# Cut 0 at the end:
					assert(sat_assignment[-1] == '0')
					sat_assignment = sat_assignment[:-1]
					break
	return sat_assignment

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
	res_ls = []
	for i in range(ls_order):
		row = []
		for j in range(ls_order):
			row.append(ls[i*ls_order + j])
		res_ls.append(row)
	assert(len(res_ls) == ls_order)
	return res_ls

def transversals_from_sat(trv_sat_assignment : list, ls_order : int, transversals_num : int):
	assert(len(trv_sat_assignment) > 0)
	assert(ls_order > 0)
	assert(transversals_num > 0)
	res_trvs = []
	trv = []
	cell_literals = []
	for lit in trv_sat_assignment:
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
			trv.append(val)
			cell_literals = []
			if len(trv) == ls_order:
				res_trvs.append(trv)
				trv = []
	assert(len(res_trvs) == transversals_num)
	return res_trvs

def are_arrays_unique(trvs : list):
    unique_trvs_set = set()
    for trv in trvs:
        unique_trvs_set.add(tuple(trv))
    return len(unique_trvs_set) == len(trvs)

def is_eo_array_correct(array : list):
	new_set = set(array)
	if len(new_set) != len(array):
		return False
	return True

### Main function:
if len(sys.argv) < 3:
	print('Usage : ls-order transversals-num [--diag] [--split]')
	print('  ls-order : order of Latin squares')
	print('  transversals-num : number of transversals to be found')
	print('  --diag   : if given, then both Latin sqaures are diagonal')
	exit(1)

print('Script ' + script + ' of version ' + version + ' is running')
if sys.argv[1] == '-v':
	exit(1)

# Parse input parameters:
ls_order = int(sys.argv[1])
transversals_num = int(sys.argv[2])
assert(ls_order > 0)
assert(transversals_num > 0)
is_diag = False
if len(sys.argv) > 2 and '--diag' in sys.argv[3:]:
	is_diag = True
is_split = False
if len(sys.argv) > 2 and '--split' in sys.argv[3:]:
	is_split = True
print('ls_order         : ' + str(ls_order))
print('transversals_num : ' + str(transversals_num))
print('is_diag          : ' + str(is_diag))
print('is_split         : ' + str(is_split))

vars_num = 0

# The first ls_order^3 variables in the CNF encode Latin square X:
ls_variables = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
# Variables for Latin squares X:
ls_vars_num = 0
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			ls_vars_num += 1
			ls_variables[i][j][k] = vars_num

clauses = []
# Latin square constraints for a Latin square:
clauses += latin_square_clauses(ls_variables, is_diag)

ls_clauses_num = len(clauses)
print(str(vars_num) + ' variables encode Latin square')
print(str(ls_clauses_num) + ' clauses encode Latin squares constraints')

if is_split:
	split_clauses = known_main_diagonal_clauses(ls_variables)
	clauses += split_clauses
	print(str(len(split_clauses)) + ' clauses encode that 2 main diagonals are known')
else:
	first_row_clauses = first_row_asending_clauses(ls_variables)
	clauses += first_row_clauses
	print(str(len(first_row_clauses)) + ' clauses encode that the first row is in the ascending order')

# More ls_order^2 variables encode each transversal.
transversals_variables = [[[0 for k in range(ls_order)] for j in range(ls_order)] for trv_index in range(transversals_num)]
trv_vars_num = 0
for trv_index in range(transversals_num):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			trv_vars_num += 1
			transversals_variables[trv_index][j][k] = vars_num
	
assert(trv_vars_num == transversals_num*ls_order*ls_order)
assert(vars_num == ls_order*ls_order*ls_order + transversals_num*ls_order*ls_order)
print(str(trv_vars_num) + ' variables encode transversals')

cpindex_variables = [[[0 for k in range(ls_order)] for j in range(ls_order)] for trv_index in range(transversals_num)]
cpindex_vars_num = 0
for trv_index in range(transversals_num):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			cpindex_vars_num += 1
			cpindex_variables[trv_index][j][k] = vars_num

assert(cpindex_vars_num == trv_vars_num)
print(str(cpindex_vars_num) + ' CP-index variables')

# Transversals constraints:
trv_clauses = transversals_clauses(ls_order, ls_variables, cpindex_variables, transversals_variables, is_diag)
print(str(len(trv_clauses)) + ' clauses encode transverals')
clauses += trv_clauses

elements = [trv_index for trv_index in range(transversals_num)]
k = 2 
# Generate all combinations of 'k' elements from 'elements'
combinations = itertools.combinations(elements, k)

combinations_num = 0
xor_vars_num = 0
xor_clauses = []
combinations_num = math.comb(len(elements), k)
print(str(combinations_num) + ' combinations of pairs of transversals')

processed_combo_num = 0
for combo in combinations:
    if processed_combo_num > 0 and processed_combo_num % 10000 == 0:
        print(str(processed_combo_num) + ' combinations processed')
    #print(combo)
    assert(len(combo) == 2)
	# Inrtoduce N^2 variables for XORs of a pair of transversals:
    xor_variables = [[0 for k in range(ls_order)] for j in range(ls_order)]
    for j in range(ls_order):
        for k in range(ls_order):
            vars_num += 1
            xor_vars_num += 1
            xor_variables[j][k] = vars_num
    combo_xor_clauses = different_transversals_clauses(ls_order, transversals_variables[combo[0]], transversals_variables[combo[1]], xor_variables)
    xor_clauses += combo_xor_clauses
    processed_combo_num += 1

assert(processed_combo_num == combinations_num)

print(str(xor_vars_num) + ' variables encode difference of transversals')
print(str(len(xor_clauses)) + ' clauses encode difference of transversals')
clauses += xor_clauses

print(str(vars_num) + ' variables in total')
print(str(len(clauses)) + ' clauses in total')

cnf_name = "transversals_n=" + str(ls_order) + '_' + str(transversals_num) + 'trvs' 
solver = 'kissat4.0.3'
res_solver_fname = 'out_n=' + str(ls_order) + '_' + str(transversals_num) + 'trvs'
if is_diag:
	cnf_name += '_diag'
	res_solver_fname += '_diag'
if is_split:
	cnf_name += '_split'
	res_solver_fname += '_split'
cnf_name += '.cnf'

print('Writing to file ' + cnf_name + ' ...')
with open(cnf_name, 'w') as ofile:
	ofile.write('p cnf ' + str(vars_num) + ' ' + str(len(clauses)) + '\n')
	for cla in clauses:
		s = ''
		for lit in cla:
			s += str(lit) + ' '
		ofile.write(s + '0\n')

command = solver + ' ' + cnf_name + ' > ' + res_solver_fname
print(command)
os.system(command)

sat_assignment = parse_sat_assignments(res_solver_fname, ls_order)
if len(sat_assignment) == 0:
	print('UNSAT')
else:
	ls_sat_assignment = sat_assignment[:ls_vars_num]
	ls = ls_from_sat(ls_sat_assignment, ls_order)
	for row in ls:
		print(row)
	trv_sat_assignment = sat_assignment[ls_vars_num:ls_vars_num + trv_vars_num]
	assert(len(trv_sat_assignment) == transversals_num*ls_order*ls_order)
	trvs_indices = transversals_from_sat(trv_sat_assignment, ls_order, transversals_num)
	assert(are_arrays_unique(trvs_indices))
	print('')
	for trv_indices in trvs_indices:
		assert(is_eo_array_correct(trv_indices))
		print(trv_indices)
	print('')
	print('Transversals :')
	trvs = []
	for trv_indices in trvs_indices:
		s = ''
		trv = []
		for i in range(len(trv_indices)):
			trv.append(ls[i][trv_indices[i]])
			s += str(ls[i][trv_indices[i]]) + ' '
		assert(is_eo_array_correct(trv))
		trvs.append(trv)
		print(s)
	assert(are_arrays_unique(trvs))
