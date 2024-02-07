# Created on: 15 Sep 2023
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# For a given CNF (that encodes the search for two orthogonal Latin squares),
# a set of cells mapping schemas (CMS), and a set of partial Latin squares, 
# for the Cartesian product of CMSs and partial Latin squares, add to the CNF
# the CMS constraints and substitute the partial Latin square.
#==============================================================================

import sys
import math

script = "cms_and_partialls_to_2mols.py"
version = "0.1.0"

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

# Generate clauses that encode the CMS condition:
def gen_cms_clauses(ls_order : int, cms : list):
  clauses_cms = []
  for i in range(ls_order):
    for j in range(ls_order):
      first_ls_cell_vars = [cnf_var_num(ls_order, 0, i, j, k) for k in range(ls_order)]
      i2 = int(cms[i][j] / ls_order)
      j2 = cms[i][j] % ls_order
      second_ls_cell_vars = [cnf_var_num(ls_order, 1, i2, j2, k) for k in range(ls_order)]
      for k in range(len(first_ls_cell_vars)):
        clauses_cms.append(str(first_ls_cell_vars[k]) + ' -' + str(second_ls_cell_vars[k]) + ' 0\n')
        clauses_cms.append('-' + str(first_ls_cell_vars[k]) + ' ' + str(second_ls_cell_vars[k]) + ' 0\n')
  assert(len(clauses_cms) == pow(ls_order,3)*2)
  return clauses_cms

print(script + ' of version ' + version + ' is running')
if len(sys.argv) == 2 and sys.argv[1] == '-v':
	exit(1)

if len(sys.argv) < 5:
	print('Usage : ls-order cnf cms-file patial-ls-file')
	print('  ls-order        : Latin square order')
	print('  cnf             : a CNF that encodes MOLS for 2 Latin squares')
	print('  cms-file        : a file with a CMS')
	print('  partial-ls-file : file with partil Latin squares')
	exit(1)

ls_order = int(sys.argv[1])
cnf_file_name = sys.argv[2]
cms_file_name = sys.argv[3]
partial_ls_file_name = sys.argv[4]

assert(ls_order > 0)
assert(cnf_file_name != '')
assert(cms_file_name != '')
assert(partial_ls_file_name != '')

print('ls_order             : ' + str(ls_order))
print('cnf_file_name        : ' + cnf_file_name)
print('cms_file_name        : ' + cms_file_name)
print('partial_ls_file_name : ' + partial_ls_file_name)

var_num = 0
clauses = []

# Read CMSs from the file:
cms_lst = []
with open(cms_file_name, 'r') as cmsf:
	lines = cmsf.read().splitlines()
	cms = []
	for line in lines:
		if 'Loading' in line or 'ESODLS' in line or len(line.split()) < 3:
			continue
		row = [int(i) for i in line.split()]
		assert(len(row) == ls_order)
		cms.append(row)
		if len(cms) == ls_order:
			cms_lst.append(cms)
			cms = []

print(str(len(cms_lst)) + ' cms were read')
print('The first CMS:')
for row in cms_lst[0]:
     print(row)
print('')
print('The last CMS:')
for row in cms_lst[-1]:
     print(row)
print('')

# Read and parse all partial Latin squares:
partial_squares = []
with open(partial_ls_file_name, 'r') as f:
    lines = f.read().splitlines()
    square = []    
    for line in lines:
      if len(line) < 2 or 'filling' in line:
        continue
      row = line.split()
      square.append(row)
      if len(square) == ls_order:
        partial_squares.append(square)
        square = []
assert(len(partial_squares) > 0)
print(str(len(partial_squares)) + ' partial LSs were read')
print('The first partial LS:')
for row in partial_squares[0]:
  print(row)
print('')
print('The last partial LS:')
for row in partial_squares[-1]:
  print(row)
print('')

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
cms_num = 0
cnf_num = 0
cla_num = 0
max_digits_cms_num = int(math.log10(len(cms_lst))) + 1
max_digits_pls_num = int(math.log10(len(partial_squares))) + 1
for cms in cms_lst:
  clauses_cms = gen_cms_clauses(ls_order, cms)
  pls_num = 0
  for pls in partial_squares:
    clauses_partial_ls = gen_partial_ls_clauses(pls)
    cms_num_str = str(cms_num)
    if len(cms_num_str) < max_digits_cms_num:
      cms_num_str = '0'*(max_digits_cms_num - len(cms_num_str)) + cms_num_str
    pls_num_str = str(pls_num)
    if len(pls_num_str) < max_digits_pls_num:
      max_digits_pls_num = '0'*(max_digits_pls_num - len(pls_num_str)) + pls_num_str
    new_cnf_name = cnf_file_name.replace('.cnf','') + '_cms=' + cms_num_str +\
      '_filling=' + pls_num_str + '.cnf'
    with open(new_cnf_name, 'w') as cnf:
      cla_num = len(clauses) + len(clauses_cms) + len(clauses_partial_ls)
      cnf.write('p cnf ' + str(var_num) + ' ' + str(cla_num) + '\n')
      for clause in clauses:
        cnf.write(clause)
      for clause in clauses_cms:
        cnf.write(clause)
      for clause in clauses_partial_ls:
        cnf.write(clause)
    pls_num += 1
  cnf_num += pls_num 
  cms_num += 1
assert(cms_num > 0)
assert(cnf_num > 0)
assert(pls_num > 0)
print(str(cnf_num) + ' CNFs were generated')
print('var num     : ' + str(var_num))
print('clauses num : ' + str(cla_num))
