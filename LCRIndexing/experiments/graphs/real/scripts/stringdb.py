import csv, sys
from itertools import count
from collections import defaultdict

fname = sys.argv[1]
edgeFname = sys.argv[2]

#with open('9606.protein.actions.v10.txt', 'r') as f:

c = count(1)
d = defaultdict(lambda: next(c))

m = ['activation', 'catalysis', 'expression', 'inhibition', 'reaction', 'binding', 'ptmod']
m = dict(zip(m, range(len(m))))

with open(fname, 'r') as f:
    with open(edgeFname, 'w') as out:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            a, b, mode, _, forward, _ = row
            if forward == 'f':
                a, b = b, a
            elif forward == 't':
                pass
            else:
                print('forward is', forward, ', unexpected')
                exit(1)

            a = d[a]
            b = d[b]
            l = m[mode]

            out.write(f'{a} {b} {l}\n')
