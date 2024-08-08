# Created on: 8 Aug 2023
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# For a given CNF (that encodes the search for MOLS),and a partial Latin
# square, add the partial Latin square constraint to the first Latin square in
# the CNF.
#==============================================================================

import sys
import math

script = "partialls_to_cnf.py"
version = "0.0.1"

# Given information on a LS's cell, return its variable in a CNF:
def cnf_var_num(ls_order : int, ls_index : int, row_index : int, \
                col_index : int, cell_val : int):
	return ls_index*pow(ls_order,3) + row_index*pow(ls_order,2) + \
    col_index*ls_order + cell_val + 1

# Generate clauses that encode a partial Latin square
def gen_partial_ls_clauses(pls : list):
  clauses_partial_ls = []
  for i in range(len(pls)):
    for j in range(len(pls[i])):
      if pls[i][j] == '-':
        continue
      var = cnf_var_num(ls_order, 0, i, j, int(pls[i][j]))
      clauses_partial_ls.append(str(var) + ' 0\n')
  assert(len(clauses_partial_ls) > 0)
  return clauses_partial_ls

print(script + ' of version ' + version + ' is running')
if len(sys.argv) == 2 and sys.argv[1] == '-v':
	exit(1)

if len(sys.argv) < 3:
	print('Usage : ls-order cnf  patial-ls-file')
	print('  ls-order        : Latin square order')
	print('  cnf             : a CNF that encodes MOLS for 2 Latin squares')
	print('  partial-ls-file : file with partil Latin squares')
	exit(1)

ls_order = int(sys.argv[1])
cnf_file_name = sys.argv[2]
partial_ls_file_name = sys.argv[3]

assert(ls_order > 0)
assert(cnf_file_name != '')
assert(partial_ls_file_name != '')

print('ls_order             : ' + str(ls_order))
print('cnf_file_name        : ' + cnf_file_name)
print('partial_ls_file_name : ' + partial_ls_file_name)

var_num = 0
clauses = []

# Read and parse a partial Latin square:
partial_squares = []
with open(partial_ls_file_name, 'r') as f:
    lines = f.read().splitlines() 
    for line in lines:
      if len(line) < 2 or 'filling' in line:
        continue
      row = line.split()
      partial_squares.append(row)
assert(len(partial_squares) == ls_order)
print('The partial LS:')
for row in partial_squares:
  print(row)

# Read clauses from a given CNF:
with open(cnf_file_name, 'r') as cnf:
	lines = cnf.readlines()
	for line in lines:
		if len(line) < 2 or line[0] == 'c':
			continue
		if line[0] == 'p':
			var_num = int(line.split(' ')[2])
		else:
			clauses.append(line)
print('clauses num   : ' + str(int(len(clauses))))
print('variables num : ' + str(var_num))
print('')
assert(var_num > 0)
assert(len(clauses) > 0)

# Generate CNFs with added CMS- and partial-clauses:

clauses_partial_ls = gen_partial_ls_clauses(partial_squares)
new_cnf_name = cnf_file_name.replace('.cnf','') + '_partialls.cnf'
print('Writing to CNF ' + new_cnf_name)
with open(new_cnf_name, 'w') as cnf:
  cla_num = len(clauses) + len(clauses_partial_ls)
  cnf.write('p cnf ' + str(var_num) + ' ' + str(cla_num) + '\n')
  for clause in clauses_partial_ls:
    cnf.write(clause)
  for clause in clauses:
    cnf.write(clause)

print('new var num     : ' + str(var_num))
print('new clauses num : ' + str(cla_num))
