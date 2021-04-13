#!/usr/bin/python3

# Usage: stdin edge file -> stdout with labels

import sys, random, math, tempfile, os, stdin
from operator import itemgetter

L = 8
for line in stdin:
    u, v = line.split()
    label = random.expovariate( 1.0 / (L / 1.7 ) )
    label = max(0,label);
    label = min(L-1,label);
    label = int(math.floor(label));
    print(f'{u} {v} {label}')
