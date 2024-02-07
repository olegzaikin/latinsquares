# Created on: 1 Sep 2023
# Author: Oleg Zaikin
# E-mail: zaikin.icc@gmail.com
#
# For a given Latin square and CMS, produce the Latin square that
# matches the CMS
#=============================================================================

import sys
import math

scriptname = 'test_ls-cms.py'
version = "0.0.1"

def print_ls(ls):
   for row in ls:
      s = ''
      for val in row[:-1]:
         s += str(val) + ' '
      s += str(row[-1])
      print(s)

if len(sys.argv) < 3:
    print('Usage: ' + scriptname + ' ls cms')
    exit(1)

print('Running ' + scriptname + ' of version ' + version)

lsname = sys.argv[1]
cmsname = sys.argv[2]
print('lsname  : ' + lsname)
print('cmsname : ' + cmsname)

ls = []
with open(lsname, 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
      ls.append([int(x) for x in line.split()])
order = len(ls)
print('LS :')
print_ls(ls)

cms = []
with open(cmsname, 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
      cms.append([int(x) for x in line.split()])
print('CMS :')
for row in cms:
    print(row)

second_ls = [[0]*order for i in range(order)]
for i in range(len(cms)):
  for j in range(len(cms[i])):
    val = cms[i][j]
    new_i = math.floor(val / order)
    new_j = val % order
    second_ls[new_i][new_j] = ls[i][j]
print('2nd LS:')
print_ls(second_ls)

for i in range(len(ls)):
   s = ''
   for j in range(len(ls[i])):
      s += str(ls[i][j]) + str(second_ls[i][j]) + ' '
   print(s)
