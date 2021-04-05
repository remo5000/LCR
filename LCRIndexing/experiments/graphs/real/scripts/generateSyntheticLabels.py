#!/usr/bin/python3

import sys, random, math, tempfile, os
from operator import itemgetter


# Usage: stdin edge file -> stdout with labels

# if(len(sys.argv) != 2):
#     print("usage python2.7 " + sys.argv[0] + " <edge_file>");
#     sys.exit(1);
# try:
#     edge_file = sys.argv[1];
# except:
#     print("Provide a filename")
#     sys.exit(1);

from sys import stdin
L = 8
for line in stdin:
    u, v = line.split()
    label = random.expovariate( 1.0 / (L / 1.7 ) )
    label = max(0,label);
    label = min(L-1,label);
    label = int(math.floor(label));
    print(f'{u} {v} {label}')


# with open(edge_file + '.tmp', 'w')  as tmp:
#     with open(edge_file, 'r') as inp:
#         for line in inp:
#                 u, v = line.split()
# 
#                 label = random.expovariate( 1.0 / (L / 1.7 ) )
#                 label = max(0,label);
#                 label = min(L-1,label);
#                 label = int(math.floor(label));
# 
#                 tmp.write(f'{u} {v} {label}\n')
# 
# with open(edge_file, 'w') as out:
#     with open(edge_file + '.tmp', 'r')  as tmp:
#         for line in tmp:
#             out.write(line)
# 
# os.remove(edge_file + '.tmp')
