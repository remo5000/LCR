#!/usr/bin/python3

import csv, sys
from itertools import count
from collections import defaultdict

from sys import stdin


c = count(1)
d = defaultdict(lambda: next(c))
m = ['activation', 'catalysis', 'expression', 'inhibition', 'reaction', 'binding', 'ptmod']
m = dict(zip(m, range(len(m))))

for line in stdin:
	if len(line.split()) == 6:
		continue
	a, b, mode, a_is_acting, _ = line.split()
	print(a, b, mode, a_is_acting)
	if a_is_acting == 't':
		pass
	elif a_is_acting == 'f':
		a, b = b, a
	else:
		print('a_is_acting is', a_is_acting, ', unexpected')
		exit(1)
	print('->', a, b, mode, a_is_acting)

	a = d[a]
	b = d[b]
	l = m[mode]

	print('->', a, b, l)

# with open(fname, 'r') as f:
#     with open(edgeFname, 'w') as out:
#         reader = csv.reader(f, delimiter='\t')
#         next(reader)
#         for row in reader:
#             a, b, mode, _, forward, _ = row
#             if forward == 'f':
#                 a, b = b, a
#             elif forward == 't':
#                 pass
#             else:
#                 print('forward is', forward, ', unexpected')
#                 exit(1)
# 
#             a = d[a]
#             b = d[b]
#             l = m[mode]
# 
#             out.write(f'{a} {b} {l}\n')
