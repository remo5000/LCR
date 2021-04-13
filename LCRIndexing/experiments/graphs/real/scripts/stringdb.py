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

	if a_is_acting == 't':
		pass
	elif a_is_acting == 'f':
		a, b = b, a
	else:
		print('a_is_acting is', a_is_acting, ', unexpected')
		exit(1)

	a = d[a]
	b = d[b]
	l = m[mode]

	print('->', a, b, l)
