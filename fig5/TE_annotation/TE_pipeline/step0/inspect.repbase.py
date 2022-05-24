def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.strip('\n').strip('>')
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	return(lines)


L = []

fadict = read_fasta('TREP.RM.combined.final')
for header, seq in fadict.items():
	bn = header.split()[0]
	L.append(bn.split("#")[1])

with open("inspect.repbase.names", 'w') as outF:
	for e in list(set(L)):
		outF.write('{}\n'.format(e))


