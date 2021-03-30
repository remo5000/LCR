m = {
'1': 0,
'.8': 1,
'.6': 2,
}

with open('out.advogato', 'r') as inp:
	with open('advogato.edge', 'w') as out:
		for line in inp:
			if '%' in line: 
				continue
			u,v,l = line.split()
			out.write(f'{u} {v} {m[l]}\n')
