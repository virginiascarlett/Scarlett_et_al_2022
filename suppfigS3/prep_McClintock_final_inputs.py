"""
McClintock pipline complained about our consensus TE library:
"Problematic symbols: ; & ( ) | * ? [ ] ~ { } < ! ^ " ' $ / + - #"

Run this script after you've already run prep_McClintock_raw_gff.py.
This script will create a new TE library and a new gff with McClintock-acceptable 
headers (no special characters).
This script will also produce a TSV file with family classifications.

run in the same directory where the consensus TE library is located, 
and specify the genome and then the library filename like so e.g.
python3 prep_McClintock_final_inputs.py ABR113 masterlib.nr.fa
"""

import re
import collections
import sys

def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.strip('\n')#.strip('>')
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	return(lines)



genome = sys.argv[1]
libname = sys.argv[2]
TEdata = read_fasta(libname)
TEnames = list(TEdata.keys())
TEnames_stripped = [ re.sub("-", "DASH", x.split('#')[0]) for x in TEnames ] #replace dashes to make consensus names unique, replace with letters
#because TE names in the output files will be messed up otherwise
TEnames_stripped = [ re.sub("_", "UNDERSCORE", x.split('#')[0]) for x in TEnames_stripped ] 
if len(TEnames_stripped) != len(set(TEnames_stripped)):
	print("{} fasta headers are not unique".format(str(len(TEnames_stripped)-len(set(TEnames_stripped)))))
	print([item for item, count in collections.Counter(TEnames_stripped).items() if count > 1])
else:
	with open('{}.clean.fa'.format(libname[:-3]), 'w') as outF:
		for n in range(len(TEnames)):
			outF.write('{}\n'.format(TEnames_stripped[n]))
			outF.write('{}\n'.format(TEdata[TEnames[n]]))

TEclassifications = {}
gfflines = []
TEcopynames = {} 
gffname = '{}_RepeatMasker_{}.gff'.format(genome, libname.split('.')[0])
with open(gffname, 'r') as mygff:
	for line in mygff.readlines():
		linedata = line.split('\t')
		familyname = re.sub("-", "DASH", linedata[8])[3:].strip('\n') #changes e.g. ID=BdisAB.trep00678\n to BdisAB.trep00678
		familyname = re.sub("_", "UNDERSCORE", familyname)
		uniquename = "{}chr{}start{}".format(familyname, linedata[0], linedata[3]) #Use TEname_chrom_start to make IDs unique
		#replace dashes so that individual TE base names match consensus names
		linedata[8] = "ID={}".format(uniquename)
		if not uniquename in TEcopynames: # If two TEs start at the same base on the same chrom, just randomly skip one of them
			TEcopynames[uniquename] = ''
			gfflines.append("\t".join(linedata))
			TEclassifications[uniquename] = familyname
		

with open('{}.clean.gff'.format(gffname[:-4]), 'w') as outgff:
	for line in gfflines:
		outgff.write("{}\n".format(line))

with open('{}_{}_McClintockfams.tsv'.format(genome, libname.split('.')[0]), 'w') as outtsv:
	for k, v in TEclassifications.items():
		outtsv.write('{}\t{}\n'.format(k, v))
