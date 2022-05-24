import timeForTE
import re


ints = []
LTRs = []
singletons = []
with open('{}.dict'.format(timeForTE.genome), 'r') as inF:
	for line in inF.readlines():
		L = line.split('\t')
		if len(L) == 2:
			ints.append(L[0])
			LTRs.append(L[1].strip('\n'))
		else:
			singletons.append(L[0].strip('\n'))

intregex = ["\-int$", "_int$", "\-I$", "_I$", "\-INT$", "_INT$"]
for myTE in singletons: #e.g. ABR113_LTRharvest00001-LTR, TREP_Ivana.trep00237-int, Repbase_IRRE_IB
	for p in intregex:
		if re.search(p, myTE):  #e.g. TREP_Ivana.trep00237-int
			stripped = re.sub(p, "", myTE) # e.g. TREP_Ivana.trep00237
			occurrences = [TE for TE in singletons if TE.startswith(stripped)]  #e.g. ['TREP_Ivana.trep00237-int', 'TREP_Ivana.trep00237-LTR']
			occurrences.remove(myTE)  #e.g. ['TREP_Ivana.trep00237-LTR']
			if len(occurrences) > 0:
				if occurrences[0].endswith('-LTR'): #double-check that we have a real int/LTR pair
					LTRs.append(occurrences[0])
					ints.append(myTE)
					break #breaks the inner for loop, i.e. 'for p in intregex:'


singletons = [TE for TE in singletons if not TE in ints]
singletons = [TE for TE in singletons if not TE in LTRs]


with open('{}.dict'.format(timeForTE.genome), 'w') as outF: #OMG overwriting the input file!
	for n in range(len(ints)):
		outF.write('{}\t{}\n'.format(ints[n].strip('\n'), LTRs[n].strip('\n')))
	for TE in singletons:
		outF.write('{}\n'.format(TE))

