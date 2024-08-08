# Created on: 8 Aug 2024
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# For a given order n, produce a CNF that encodes the search for 3 MOLS
# (mutually orthogonal Latin squares) of order n. The SAT encoding is CP-index
# proposed in
#   Noah Rubin, Curtis Bright, Brett Stevens, Kevin Cheung. Integer and 
#   Constraint Programming Revisited for Mutually Orthogonal Latin Squares 
#   // In AAAI 2022.
# For the mentioned paper, the CP-index SAT encoding was implemented in C++ in
#   https://github.com/noahrubin333/CP-IP 
#=============================================================================

import sys

script = "cp-index_sat_3mols.py"
version = "0.0.1"

# Return clauses that encode the AtMostOne constraint:
def at_most_one_clauses(vars : list):
	assert(len(vars) > 1)
	assert(len(set(vars)) == len(vars))
	res_clauses = []
	for i in range(len(vars)):
		for j in range(i+1, len(vars)):
			res_clauses.append([-vars[i], -vars[j]])
	return res_clauses

# Get clauses encoding the ExactlyOne constraint for a set of variables:
def exactly_one_clauses(vars : list):
	assert(len(vars) > 1)
	assert(len(set(vars)) == len(vars))
	res_clauses = at_most_one_clauses(vars)
	res_clauses.append(vars) # AtLeastOne constraint
	return res_clauses

# Latin square constraints:
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

# Orthogonality is encoded via element indexing constraint Z_i[X_ij] = Yij.
# This constraint for a certain value in a Z's cell is that
# if cell X[i][j] has value k and Y[i][j] has value l,
# then in i-th row of Z value l is in column k.
# In other words, 
# (X[i][j][k] & Y[i][j][l]) -> Z[i][k][l]
# (Y[i][j][l] & Z[i][k][l]) -> X[i][j][k]
# (X[i][j][k] & Z[i][k][l]) -> Y[i][j][l]
# This is encoded via three clauses:
# [-Y[i][j][l], -X[i][j][k], Z[i][k][l]]
# [-Y[i][j][l], -Z[i][k][l], X[i][j][k]]
# [-Z[i][k][l], -X[i][j][k], Y[i][j][l]]
# The AllDifferent constraint for Z's columns encures the orhtogonality,
# yet the AllDifferent for Z's rows redundant but may help a solver.
def orthogonality_clauses(LS1 : list, LS2 : list, Orth_LS : list):
	res_clauses = []
	ls_order = len(LS1)
	assert(ls_order > 0 and ls_order == len(LS2) and ls_order == len(Orth_LS))
	for i in range(ls_order):
		for j in range(ls_order):
			for k in range(ls_order):
				for l in range(ls_order):
					res_clauses.append([-LS1[i][j][k], -LS2[i][j][l], Orth_LS[i][k][l]])
					res_clauses.append([-LS2[i][j][l], -Orth_LS[i][k][l], LS1[i][j][k]])
					res_clauses.append([-LS1[i][j][k], -Orth_LS[i][k][l], LS2[i][j][l]])
	return res_clauses

### Main function:
if len(sys.argv) < 3:
	print('Usage : ls-order cnf-name [--diag]')
	print('  ls-order : order of Latin squares')
	print('  cnf-name : name of output CNF')
	print('  --diag   : if given, then both Latin sqaures are diagonal')
	exit(1)

print('Script ' + script + ' of version ' + version + ' is running')
if sys.argv[1] == '-v':
	exit(1)

# Parse input parameters:
ls_order = int(sys.argv[1])
cnf_name = sys.argv[2]
is_diag = False
if len(sys.argv) > 3 and '--diag' in sys.argv[3:]:
	is_diag = True
print('ls_order : ' + str(ls_order))
print('cnf_name : ' + cnf_name) 
print('is_diag  : ' + str(is_diag))

vars_num = 0

# The first ls_order^3 variables in the CNF encode Latin square X:
X = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
# Variables for Latin squares X:
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			X[i][j][k] = vars_num
			
# Next ls_order^3 variables encode Latin square Y:
Y = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			Y[i][j][k] = vars_num

# Next ls_order^3 variables encode Latin square Z:
Z = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			Z[i][j][k] = vars_num

# ls_order^3 for a Latin square that ensures orthogonality for X and Y:
OrthLS_XY = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			OrthLS_XY[i][j][k] = vars_num

# ls_order^3 for a Latin square that ensures orthogonality for X and Z:
OrthLS_XZ = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			OrthLS_XZ[i][j][k] = vars_num

# ls_order^3 for a Latin square that ensures orthogonality for Y and Z:
OrthLS_YZ = [[[0 for k in range(ls_order)] for j in range(ls_order)] for i in range(ls_order)]
for i in range(ls_order):
	for j in range(ls_order):
		for k in range(ls_order):
			vars_num += 1
			OrthLS_YZ[i][j][k] = vars_num

clauses = []
# Latin square constraints for all Latin squares:
clauses += latin_square_clauses(X, is_diag)
clauses += latin_square_clauses(Y, is_diag)
clauses += latin_square_clauses(Z, is_diag)
# Diagonal constraints are not needed for orthogonal-Latin-squares:
clauses += latin_square_clauses(OrthLS_XY, False)
clauses += latin_square_clauses(OrthLS_XZ, False)
clauses += latin_square_clauses(OrthLS_YZ, False)

ls_clauses_num = len(clauses)
print(str(ls_clauses_num) + ' clauses encode Latin squares constraints')

# Orthogonality constraints:
clauses += orthogonality_clauses(X, Y, OrthLS_XY)
clauses += orthogonality_clauses(X, Z, OrthLS_XZ)
clauses += orthogonality_clauses(Y, Z, OrthLS_YZ)

print(str(len(clauses) - ls_clauses_num) + ' clauses encode orthogonality constraints')
print(str(len(clauses)) + ' clauses in total')

print('Writing to file ' + cnf_name + ' ...')
with open(cnf_name, 'w') as ofile:
	ofile.write('p cnf ' + str(vars_num) + ' ' + str(len(clauses)) + '\n')
	for cla in clauses:
		s = ''
		for lit in cla:
			s += str(lit) + ' '
		ofile.write(s + '0\n')
