"""
A janky script to remove simple repeats from a RepeatMasker database
based on the phrases 'Simple_repeat' or 'Low_complexity' appearing in 
the fasta header. 
"""

fobj = open('TREP.RM.combined', 'r')
alldata = {}
current_elem = ''
current_elem_data = []
for line in fobj.readlines():
	if line.startswith('>'):
		if current_elem!='':
			alldata[current_elem] = ''.join(current_elem_data)
		current_elem = line
		current_elem_data = []
	elif not line.startswith('>'):
		current_elem_data.append(line)


alldata[current_elem] = ''.join(current_elem_data)

fobj.close()

finaldata={}
for k, v in alldata.items():
	if not 'Simple_repeat' in k:
		finaldata[k] = v
	
finalfinaldata={}
for k, v in finaldata.items():
	if not 'Low_complexity' in k:
		finalfinaldata[k] = v

outF = open('TREP.RM.combined', 'w')
for k, v in finalfinaldata.items():
	outF.write(k)
	outF.write(v)

outF.close()

