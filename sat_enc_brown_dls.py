# Created on: 1 Nov 2024
# Author: Oleg Zaikin
# E-mail: oleg.zaikin@icc.ru
#
# SAT encoding for enumeration of all Brown squares of given order.
# See the description here https://oeis.org/A339305
#==============================================================================

import sys
import itertools
import os

script = "sat_enc_brown_dls.py"
version = "0.0.3"

def print_row(arr : list):
    assert(len(arr) > 0)
    s = ''
    for x in arr[:-1]:
        s += str(x) + " "
    s += str(arr[-1])
    print(s)

def clause_to_str(clause : list):
    s = ''
    for lit in clause:
        s += str(lit) + ' '
    return s + '0'

# Clauses for the AtMostOne constraint:
def at_most_one_clauses(vars : list):
	assert(len(vars) > 1)
	assert(len(set(vars)) == len(vars))
	res_clauses = []
	for i in range(len(vars)-1):
		for j in range(i+1, len(vars)):
			res_clauses.append([-vars[i], -vars[j]])
	return res_clauses

# Clauses for the ExactlyOne constraint for a set of variables:
def exactly_one_clauses(vars : list):
	assert(len(vars) > 1)
	assert(len(set(vars)) == len(vars))
	res_clauses = at_most_one_clauses(vars)
	res_clauses.append(vars) # AtLeastOne constraint
	return res_clauses

# Clauses for a Latin square:
def latin_square_clauses(T : list, is_diag : bool):
	ls_order = len(T)
	assert(ls_order > 0)
	res_clauses = []
	# Constraints on rows, columns, and values are obligatory:
	for i in range(ls_order):
		for j in range(ls_order):
			# Each square's cell contains exactly one value 0..n-1:
			res_clauses += exactly_one_clauses([T[i][j][k] for k in range(ls_order)])
			# Each value occurs exactly once in each row:
			res_clauses += exactly_one_clauses([T[i][k][j] for k in range(ls_order)])
			# Each value occurs exactly once in each column:
			res_clauses += exactly_one_clauses([T[k][i][j] for k in range(ls_order)])
	# Constraints on main diagonal and antidiagonal are optional:
	if is_diag:
		# Main diagonal:
		for i in range(ls_order):
			res_clauses += exactly_one_clauses([T[k][k][i] for k in range(ls_order)])
		# Main antidiagonal:
		for i in range(ls_order):
			res_clauses += exactly_one_clauses([T[ls_order-k-1][k][i] for k in range(ls_order)])
	return res_clauses

# Clauses and variables for diagonal Latin square:
def diag_latin_square(n : int, is_first_row_ascending : bool):
    assert(n > 0)
    # n^3 variables in the CNF encode a Latin square:
    X = [[[0 for k in range(n)] for j in range(n)] for i in range(n)]
    assert(len(X) == n)
    # Variables for Latin squares X:
    vars_num = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                vars_num += 1 # the first variable is 1
                X[i][j][k] = vars_num 
    # First row in ascending order:
    first_row_ascending_clauses = []
    if is_first_row_ascending:
        for j in range(n):
            first_row_ascending_clauses.append([X[0][j][j]])
        #print(str(len(first_row_ascending_clauses)) + ' clauses encode the first row in ascending')
    # Diagonal Latin square:
    dls_clauses = latin_square_clauses(X, True)
    assert(len(dls_clauses) > 0)
    clauses = []
    for cla in first_row_ascending_clauses:
        clauses.append(cla)
    for cla in dls_clauses:
        clauses.append(cla)
    assert(len(clauses) > 0)
    return vars_num, X, clauses

# Clauses for equality of two variables:
def equality_clauses(x : int, y : int):
    assert(x > 0 and y > 0 and x != y)
    clauses = []
    clauses.append([x, -y])
    clauses.append([-x, y])
    return clauses

# Clauses for the horizontal symmetry of a Latin square:
def horizontal_symmetry(X : list, n : int):
    assert(n > 0 and n % 2 == 0)
    assert(len(X) == n)
    horiz_sym_clauses = []
    half_n = int(n/2)
    for i in range(n): # for each row
        for j in range(half_n): # for each column 0 .. (n/2)-1
            for k in range(n): # for each cell's value
                eq_clauses = equality_clauses(X[i][j][k], X[i][n-1-j][n-1-k]) # x[i][j] == n-1-x[i][n-1-j]
                for cla in eq_clauses:
                    horiz_sym_clauses.append(cla)
    return horiz_sym_clauses

# Clauses for the vertical symmetry of a Latin square:
def vertical_symmetry(X : list, n : int, cover : list):
    assert(n > 0 and n % 2 == 0)
    assert(len(X) == n)
    assert(n > 0)
    half_n = int(n/2)
    assert(len(cover) == half_n)
    vertic_sym_clauses = []
    for i in range(half_n): # for each row 0 .. (n/2)-1
        for j in range(n): # for each column
            for inverse_columns_indices in cover:
                assert(len(inverse_columns_indices) == 2)
                eq_clauses = equality_clauses(X[i][j][inverse_columns_indices[0]], X[n-1-i][j][inverse_columns_indices[1]])
                for cla in eq_clauses:
                    vertic_sym_clauses.append(cla)
    return vertic_sym_clauses

def gen_covers(n : int):
    num_row_pairs = int(n*(n-1)/2)
    print("num_row_pairs : " + str(num_row_pairs))
    elements = [i for i in range(n)]
    pair_combinations = []
    for comb in itertools.combinations(elements, 2):
        pair_combinations.append(comb)
    assert(len(pair_combinations) == num_row_pairs)
    # Here a cover is a set of n/2 pairs which contain all n rows:
    covers = []
    comb_num = 0
    for comb in itertools.combinations(pair_combinations, int(n/2)):
        comb_num += 1
        row_indices_set = set()
        #print(comb)
        for pair in comb:
            row_indices_set.add(pair[0])
            row_indices_set.add(pair[1])
        assert(len(row_indices_set) <= n)
        if len(row_indices_set) != n:
            continue
        else:
            covers.append(comb)
    assert(len(covers) > 0)
    print(str(len(covers)) + " covers of five pairs out of " + str(comb_num) + " are possible")
    print('The first covers are :')
    print(covers[0])
    if len(covers) > 1:
        print(covers[1])
    if len(covers) > 2:
        print('The last cover is :')
        print(covers[-1])
    return covers

def generate_cnf_brown_dls_horiz_sym(n : int, vars_num : int, X : list, dls_clauses : list, cvr : list, cvr_indx : int):
    assert(n > 0 and n % 2 == 0)
    assert(vars_num > 0)
    assert(len(X) > 0)
    assert(len(dls_clauses) > 0)
    assert(len(cvr) == int(n/2))
    # Add constraints for horizontal symmetry:
    horiz_sym_clauses = horizontal_symmetry(X, n)
    assert(len(horiz_sym_clauses) > 0)
    #print(str(len(horiz_sym_clauses)) + ' clauses encode horizontal symmetry')
    # Constraints for inverse rows:
    all_eq_clauses = []
    for inverse_rows_indices in cvr:
        for j in range(n):
            for k in range(n):
                eq_clauses = equality_clauses(X[inverse_rows_indices[0]][j][k], X[inverse_rows_indices[1]][n-1-j][k])
                for cla in eq_clauses:
                    all_eq_clauses.append(cla)
    # Write to CNF:
    cnf_fname = "dls_n" + str(n) + "_first_row_ascending_horiz_sym_inverse_rows_" + str(cvr_indx) + ".cnf"
    #print('Writing a CNF for DLS with horizontal symmetry and inverse rows to the file ' + cnf_fname)
    cla_num = len(dls_clauses) + len(horiz_sym_clauses) + len(all_eq_clauses)
    with open(cnf_fname, 'w') as cnf_f:
        cnf_f.write('p cnf ' + str(vars_num) + ' ' + str(cla_num) + '\n')
        for cla in dls_clauses:
            cnf_f.write(clause_to_str(cla) + '\n')
        for cla in horiz_sym_clauses:
            cnf_f.write(clause_to_str(cla) + '\n')
        for cla in all_eq_clauses:
            cnf_f.write(clause_to_str(cla) + '\n')
    return cnf_fname

def generate_cnf_brown_dls_vertic_sym(n : int, vars_num : int, X : list, dls_clauses : list, cvr : list, cvr_indx : int):
    assert(n > 0 and n % 2 == 0)
    assert(vars_num > 0)
    assert(len(X) > 0)
    assert(len(dls_clauses) > 0)
    assert(len(cvr) == int(n/2))
    # Vartical symmetry for a given cover of pairs of column indices:
    vertic_sym_clauses = vertical_symmetry(X, n, cvr)
    #print(str(len(vertic_sym_clauses)) + ' clauses encode vertical symmetry')
    # Inverse-columns for a given cover:
    all_eq_clauses = []
    for inverse_columns_indices in cvr:
        for i in range(n):
            for k in range(n):
                eq_clauses = equality_clauses(X[i][inverse_columns_indices[0]][k], X[n-1-i][inverse_columns_indices[1]][k])
                for cla in eq_clauses:
                    all_eq_clauses.append(cla)
    cnf_fname = "dls_n" + str(n) + "_first_row_ascending_vertic_sym_inverse_columns_" + str(cvr_indx) + ".cnf"
    #print('Writing a CNF for DLS with vertical symmetry and inverse columns to the file ' + cnf_fname)
    cla_num = len(dls_clauses) + len(vertic_sym_clauses) + len(all_eq_clauses)
    with open(cnf_fname, 'w') as cnf_f:
        cnf_f.write('p cnf ' + str(vars_num) + ' ' + str(cla_num) + '\n')
        for cla in dls_clauses:
            cnf_f.write(clause_to_str(cla) + '\n')
        for cla in vertic_sym_clauses:
            cnf_f.write(clause_to_str(cla) + '\n')
        for cla in all_eq_clauses:
            cnf_f.write(clause_to_str(cla) + '\n')
    return cnf_fname

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

def ls_str_from_sat(sat_assign : list, n : int):
	assert(len(sat_assign)  == n*n*n)
	assert(n > 0 and n < 11 and n % 2 == 0)
	ls_str = ''
	cell_literals = []
	for lit in sat_assign:
		cell_literals.append(lit)
		# Literals of a current cell are collected:
		val = '-'
		neg_lit_num = 0
		if len(cell_literals) == n:
			for j in range(n):
				if int(cell_literals[j]) > 0:
					val = str(j)
				else:
					neg_lit_num += 1
			assert(neg_lit_num == n - 1)
			#assert(val >= 0 and val <= n - 1)
			ls_str += val
			cell_literals = []
	assert(len(ls_str) == n * n)
	return ls_str

def parse_latin_squares_from_sat_assign(fname : str, n : int):
    assert(n > 0 and n < 11 and n % 2 == 0)
    is_started_sol = False
    ls_lst = []
    cur_assignment = []
    with open(fname, 'r') as ifile:
        lines = ifile.read().splitlines()
		# Process divided lines with sat assignments:
        lines = merge_assign_lines(lines)
        for line in lines:
            if len(line) < 2:
                continue
            if 's SATISFIABLE' in line or 'c Answer: ' in line:
                is_started_sol = True
                cur_assignment = []
                continue
            if is_started_sol:
				# If clasp solver, no solution after 's SATISFIABLE':
                if line[0] != 'v' or line[1] != ' ':
                    continue
				# Cut 'v' at the beginning:
                line = line[1:]
				# Read literals
                literals = line.split()
                assert(ascending_order(literals) or literals[-1] == '0')
                assert(len(literals) > 0)
                for lit in literals:
                    cur_assignment.append(lit)
				# If the last line of a solution is found, stop reading:
                if literals[-1] == '0':
                    is_started_sol = False
					# Cut 0 at the end:
                    assert(cur_assignment[-1] == '0')
                    cur_assignment = cur_assignment[:-1]
					# The number of variables must be 3*N^3:
					#print(str(len(cur_assignment)))
                    ls_str = ls_str_from_sat(cur_assignment, n)
                    ls_lst.append(ls_str)
    return ls_lst

# num of unknowns rows, elements, and variants of values:
#var_num_after_propag = int((n/2 - 1) * n/2 * (n-2))
#print("var_num_after_propag :")
#print(str(var_num_after_propag))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage : " + script + " Brown-DLS-order")
        exit(1)

    n = int(sys.argv[1])
    assert(n > 1 and n < 11 and n % 2 == 0)
    print('n : ' + str(n))

    # Make a CNF that encodes a diagonal Latin square of order n:
    vars_num, X, dls_clauses = diag_latin_square(n, True)
    print(str(vars_num) + ' Boolean variables')
    print(str(len(dls_clauses)) + ' clauses encode DLS with the first row in ascending order')

    cnf_dls_fname = "dls_n" + str(n) + "_first_row_ascending.cnf"
    print('Writing a CNF for DLS to the file ' + cnf_dls_fname)
    with open(cnf_dls_fname, 'w') as cnf_dls_f:
        cnf_dls_f.write('p cnf ' + str(vars_num) + ' ' + str(len(dls_clauses)) + '\n')
        for cla in dls_clauses:
            cnf_dls_f.write(clause_to_str(cla) + '\n')

    # Generate all covers of fives of pairs: 
    covers = gen_covers(n)

    # For each cover, generate a horizontal-symmetry and a vertical-symmetry CNF:
    all_ls_set = set()
    for i in range(len(covers)):
        horiz_sym_cnf_fname = generate_cnf_brown_dls_horiz_sym(n, vars_num, X, dls_clauses, covers[i], i)
        # Find all satisfying assignments:
        out_fname = 'out_' + horiz_sym_cnf_fname.split('.cnf')[0]
        os.system('clasp --models 0 --enum-mode=bt --configuration=crafty ./' + horiz_sym_cnf_fname + ' > ' + out_fname)
        ls_lst = parse_latin_squares_from_sat_assign(out_fname, n)
        for ls in ls_lst:
            all_ls_set.add(ls)
        os.system('rm ' + horiz_sym_cnf_fname)
        os.system('rm ' + out_fname)
        #
        vertic_sym_cnf_fname = generate_cnf_brown_dls_vertic_sym(n, vars_num, X, dls_clauses, covers[i], i)
        out_fname = 'out_' + vertic_sym_cnf_fname.split('.cnf')[0]
        os.system('clasp --models 0 --enum-mode=bt --configuration=crafty ./' + vertic_sym_cnf_fname + ' > ' + out_fname)
        ls_lst = parse_latin_squares_from_sat_assign(out_fname, n)
        for ls in ls_lst:
            all_ls_set.add(ls)
        os.system('rm ' + vertic_sym_cnf_fname)
        os.system('rm ' + out_fname)
        print('processed ' + str(i+1) + ' covers out of ' + str(len(covers)))

    print(str(len(all_ls_set)) + ' Brown DLS of order ' + str(n))
    #for ls in all_ls_set:
    #    print(ls)
