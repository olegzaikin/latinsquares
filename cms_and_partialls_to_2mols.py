# Created on: 15 Sep 2023
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# For a given CNF (that encodes the search for two orthogonal Latin squares),
# a cells mapping schema (CMS), and a partial Latin square, add to the CNF
# the CMS condition and substitute the partial Latin square.
#==============================================================================

# TODO:
# 0. More than one CMS in a file

import sys

script = "cms_and_partialls_to_2mols.py"
version = "0.0.2"

# Given information on a LS's cell, return its variable in a CNF:
def cnf_var_num(ls_order : int, ls_index : int, row_index : int, \
                col_index : int, cell_val : int):
	return ls_index*pow(ls_order,3) + row_index*pow(ls_order,2) + \
    col_index*ls_order + cell_val + 1

print('Script ' + script + ' of version ' + version + ' is running')
if sys.argv[1] == '-v':
	exit(1)

if len(sys.argv) < 5:
	print('Usage : cnf cms-file patial-ls-file ls-order')
	print('  cnf             : a CNF that encodes MOLS for 2 Latin squares')
	print('  cms-file        : a file with a CMS')
	print('  partial-ls-file : file with partil Latin squares')
	print('  ls-order        : Latin square order')
	exit(1)

cnf_file_name = sys.argv[1]
cms_file_name = sys.argv[2]
partial_ls_file_name = sys.argv[3]
ls_order = int(sys.argv[4])

assert(ls_order > 0)

print('cnf_file_name        : ' + cnf_file_name)
print('cms_file_name        : ' + cms_file_name)
print('partial_ls_file_name : ' + partial_ls_file_name)
print('ls_order             : ' + str(ls_order))

var_num = 0
clauses = []

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

# Read a CMS from the file:
cms = []
with open(cms_file_name, 'r') as cmsf:
	lines = cmsf.read().splitlines()
	for line in lines:
	    row = [int(i) for i in line.split(' ')]
	    assert(len(row) == ls_order)
	    cms.append(row)
print('CMS :')
assert(len(cms) == ls_order)
for row in cms:
     print(row)
print('')

# Generate clauses that encode the CMS condition:
clauses_cms = []
for i in range(ls_order):
	for j in range(ls_order):
		first_dls_cell_vars = [cnf_var_num(ls_order, 0, i, j, k) for k in range(ls_order)]
		i2 = int(cms[i][j] / ls_order)
		j2 = cms[i][j] % ls_order
		second_dls_cell_vars = [cnf_var_num(ls_order, 1, i2, j2, k) for k in range(ls_order)]
		for k in range(len(first_dls_cell_vars)):
			clauses_cms.append(str(first_dls_cell_vars[k]) + ' -' + str(second_dls_cell_vars[k]) + ' 0\n')
			clauses_cms.append('-' + str(first_dls_cell_vars[k]) + ' ' + str(second_dls_cell_vars[k]) + ' 0\n')
assert(len(clauses_cms) == pow(ls_order,3)*2)

# Read and parse all partial Latin squares:
partial_squares = []
with open(partial_ls_file_name, 'r') as f:
    lines = f.read().splitlines()
    square = []    
    for line in lines:
      if len(line) < 2 or 'filling' in line:
        continue
      row = line.split(' ')
      if len(square) == ls_order:
        partial_squares.append(square)
        square = [row]
        #print(str(len(partial_squares)))
      else:
        square.append(row)
     # print(line)
    partial_squares.append(square)
assert(len(partial_squares) > 0)
print(str(len(partial_squares)) + ' partial LSs were read')
#for pls in partial_squares:
#  for row in pls:
#    print(row)
#  print('\n')
#exit(1)

# Generate CNFs with added CMS- and partial-clauses:
ls_num = 0
cla_num = 0
for pls in partial_squares:
  clauses_partial_ls = []
  for i in range(len(pls)):
    for j in range(len(pls[i])):
      if pls[i][j] == '-':
        continue
      var = cnf_var_num(ls_order, 0, i, j, int(pls[i][j]))
      clauses_partial_ls.append(str(var) + ' 0\n')
  #print(clauses_partial_ls)
  with open(cnf_file_name.replace('.cnf','') + '_cms=' + cms_file_name + '_filling=' + str(ls_num) + '.cnf', 'w') as cnf:
    cla_num = len(clauses) + len(clauses_partial_ls) + len(clauses_cms)
    cnf.write('p cnf ' + str(var_num) + ' ' + str(cla_num) + '\n')
    for clause in clauses_partial_ls:
      cnf.write(clause)
    for clause in clauses_cms:
      cnf.write(clause)
    for clause in clauses:
      cnf.write(clause)
  ls_num += 1
assert(ls_num > 0)
print(str(ls_num) + ' CNFs were generated')
print('clauses num : ' + str(cla_num))
