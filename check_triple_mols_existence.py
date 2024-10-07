# Created on: 4 Oct 2024
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# Given a file with Latin squares, check if a triple of MOLS can be formed
# based on them.
#=============================================================================

import sys

script = "check_triple_mols_existence.py"
version = "0.0.1"

def count_orthogonal_cells(ls1 : list, ls2 : list):
	ordered_pairs_set = set()
	for i in range(len(ls1)):
		for j in range(len(ls1[i])):
			ordered_pairs_set.add(ls1[i][j] + ls2[i][j])
	#print(ordered_pairs_set)
	return len(ordered_pairs_set)

### Main function:
print(script + ' of version ' + version + ' is running')
if len(sys.argv) == 2 and sys.argv[1] == '-v':
	exit(1)

if len(sys.argv) < 2:
	print('Usage : ' + script + ' ls-file')
	exit(1)

fname = sys.argv[1]
latin_squares = []
with open(fname, 'r') as f:
	lines = f.read().splitlines()
	for line in lines:
		line = line.rstrip()
		#print(line)
		assert(len(line) == 100)
		ls = []
		for i in range(10):
			row = []
			for j in range(10):
				row.append(line[i*10 + j])
			#print(row)
			ls.append(row)
		latin_squares.append(ls)
		#print('')
print(str(len(latin_squares)) + ' Latin squares')

ls_num = len(latin_squares)
orth_pairs_num = 0
max_orth_num = -1
for i in range(ls_num - 1):
	if i > 0 and i % 100 == 0:
		print(str(i) + ' squares processed')
	for j in range(i+1, ls_num):
		orth_num = count_orthogonal_cells(latin_squares[i], latin_squares[j])
		if orth_num > max_orth_num:
			max_orth_num = orth_num
			print('max_orth_num : ' + str(max_orth_num))
		if orth_num == 100:
			orth_pairs_num += 1
			print(str(orth_pairs_num) + ' orthogonal pairs so far.')
			#print(str(orth_num))
