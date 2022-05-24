sequences = {}
TEorder = []

with open('trep-db_nr_Rel-16.fasta', 'r') as original:
	currentHeader = ''
	for line in original.readlines():
		if line == '\n':
			pass
		else:
			if line.startswith('>'):
				currentHeader = line.strip('\n')
				sequences[currentHeader] = ''
				TEorder.append(currentHeader)
			else:
				sequences[currentHeader] += line.strip('\n')

outF = open('trep-db_nr_Rel-16.nulls-removed.fasta', 'w')
for TE in TEorder:
	if sequences[TE] == 'NULL':
		pass
	else:
		outF.write('{}\n'.format(TE))
		outF.write('{}\n'.format(sequences[TE]))
outF.close()

