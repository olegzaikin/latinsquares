# SAT encoding for enumeration of all Brown squares of given order.

import sys
import itertools

script = "sat_enc_brown_dls.py"
version = "0.0.1"

def print_row(arr : list):
    assert(len(arr) > 0)
    s = ''
    for x in arr[:-1]:
        s += str(x) + " "
    s += str(arr[-1])
    print(s)

def generate_cnf_brown_dls(n : int, cvr : tuple):
    #print(cvr)
    assert(len(cvr) > 0)
    first_row = [i for i in range(n)]
    print("first_row :")
    print_row(first_row)
    pair_of_first_row = [n-1-i for i in range(n)]
    print("pair_of_first_row :")
    print_row(pair_of_first_row)
    pair_of_first_row_index = cvr[0][1]
    print("pair_of_first_row_index : " + str(pair_of_first_row_index))
    var_num = int(n * n * n / 4)
    print("var_num :")
    print(str(var_num))
    # num of unknowns rows, elements, and variants of values:
    var_num_after_propag = int((n/2 - 1) * n/2 * (n-2))
    print("var_num_after_propag :")
    print(str(var_num_after_propag))

if len(sys.argv) < 2:
    print("Usage : " + script + " Brown-DLS-order")
    exit(1)

n = int(sys.argv[1])
assert(n > 1 and n < 11 and n % 2 == 0)
print('n : ' + str(n))

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
    for pair in comb:
        row_indices_set.add(pair[0])
        row_indices_set.add(pair[1])
    assert(len(row_indices_set) <= n)
    if len(row_indices_set) != n:
        continue
    else:
        covers.append(comb)
        #print(comb)

assert(len(covers) > 2)

print(str(len(covers)) + " covers out of " + str(comb_num) + " are possible")

print('The first 3 covers are :')
for i in range(3):
    print(covers[i])
print('The last cover is :')
print(covers[-1])

# For each cover, generate a CNF:
for cvr in covers:
    generate_cnf_brown_dls(n, cvr)
