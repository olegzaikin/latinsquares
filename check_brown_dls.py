# Created on: 1 Nov 2024
# Author: Oleg Zaikin
# E-mail: oleg.zaikin@icc.ru
#
# Check how many DLS from a given file are Brown.
#==============================================================================

import sys

script = "check_brown_dls.py"
version = "0.0.1"

def horiz_sym_row_inverse(ls : list, n : int):
    assert(len(ls) == n)
    assert(len(ls) > 0)
    assert(n > 0 and n < 11)
    half_n = int(n/2)
    for i in range(n):
        for j in range(half_n):
            if ls[i][j] != n-1-ls[i][n-1-j]:
                return False
    pair_rows_indices = [0]*n
    for i in range(n-1):
        if pair_rows_indices[i] == 1:
            continue
        for j in range(i+1, n):
            if pair_rows_indices[j] == 1:
                continue
            #print(str(i) + ' ' + str(j))
            assert(i != j)
            assert(ls[i] != ls[j])
            if ls[i] == ls[j][::-1]:
                pair_rows_indices[i] += 1
                pair_rows_indices[j] += 1
                break
    #print(pair_rows_indices)
    for x in pair_rows_indices:
        if x == 0:
            return False
    return True

def reverse_columns(ls : list, i : int, j : int):
    n = len(ls)
    assert(i != j and i < n and j < n)
    for row_index in range(n):
        if ls[row_index][i] != ls[n-1-row_index][j]:
            return False
    return True

def vertic_sym_row_inverse(ls : list, n : int):
    assert(len(ls) > 0)
    assert(n > 0 and n < 11)
    assert(n % 2 == 0)
    half_n = int(n/2)
    corresp_arr = [-1 for _ in range(n)]
    for i in range(n):
        for j in range(half_n):
            #print(str(i) + " " + str(j))
            if corresp_arr[ls[i][j]] == -1 and corresp_arr[ls[n-1-i][j]] == -1:
                corresp_arr[ls[i][j]] = ls[n-1-i][j]
                corresp_arr[ls[n-1-i][j]] = ls[i][j]
            elif (corresp_arr[ls[i][j]] == -1 and corresp_arr[ls[n-1-i][j]] != -1) or \
                 (corresp_arr[ls[i][j]] != -1 and corresp_arr[ls[n-1-i][j]] == -1):
                return False
            elif ls[n-1-i][j] != corresp_arr[ls[i][j]]:
                return False
    print_ls(ls)
    print('')
    print(corresp_arr)
    print('')
    print('')
    pair_columns_indices = [0 for _ in range(n)]
    for i in range(n-1):
        if pair_columns_indices[i] == 1:
            continue
        for j in range(i+1, n):
            if pair_columns_indices[j] == 1:
                continue
            #print(str(i) + ' ' + str(j))
            if reverse_columns(ls, i, j):
                pair_columns_indices[i] += 1
                pair_columns_indices[j] += 1
                break
    #print(pair_columns_indices)
    for x in pair_columns_indices:
        if x == 0:
            return False
    return True

def print_ls(ls : list):
    for row in ls:
        s = ''
        for x in row:
            s += str(x) + ' '
        print(s)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage : " + script + " file-with-dls dls-order")
        exit(1)

    fname = sys.argv[1]
    n = int(sys.argv[2])
    assert(n > 1 and n < 11 and n % 2 == 0)
    print('fname : ' + fname)
    print('n     : ' + str(n))

    dls_num = 0
    brown_dls_horiz_sym = 0
    brown_dls_vertic_sym = 0
    brown_dls_double_sym = 0
    with open(fname, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            #print(line)
            assert(len(line) == n*n)
            ls = []
            for i in range(n):
                row = []
                for j in range(n):
                    row.append(int(line[i*n + j]))
                ls.append(row)
            assert(len(ls) == n)
            #print_ls(ls)
            dls_num += 1
            if horiz_sym_row_inverse(ls, n):
                brown_dls_horiz_sym += 1
            if vertic_sym_row_inverse(ls, n):
                brown_dls_vertic_sym += 1
            if horiz_sym_row_inverse(ls, n) and vertic_sym_row_inverse(ls, n):
                brown_dls_double_sym += 1

    print('dls_num : ' + str(dls_num))
    #assert(brown_dls_horiz_sym + brown_dls_vertic_sym + brown_dls_double_sym == dls_num)
    print('brown_dls_horiz_sym  : ' + str(brown_dls_horiz_sym))
    print('brown_dls_vertic_sym : ' + str(brown_dls_vertic_sym))
    print('brown_dls_double_sym : ' + str(brown_dls_double_sym))
