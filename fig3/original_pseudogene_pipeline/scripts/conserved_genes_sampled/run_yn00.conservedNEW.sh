#!/bin/bash
#Run this script from two dirs above candidatePseudogeneRegions/
#need yn00.ctl in the same dir as this script
#this script will change the first two lines of it to match the current gene
for fname in $( ls candidatePseudogeneRegions/*.readyForPAML.fasta ); do 
	line="seqfile = $fname"
	#sed -i "1s/.*/$line/" yn00.ctl
	sed -i "1c\\$line" yn00.ctl
	geneID=${fname%.readyForPAML.fasta} #looks like sample1/candidatePseudogeneRegions/88
	line="outfile = $geneID.yn"
	#sed -i "2s/.*/$line/" yn00.ctl
	sed -i "2c\\$line" yn00.ctl
	/global/cfs/cdirs/plantbox/hybridum/software/paml4.9j/bin/yn00 >> yn00.stdout
done
