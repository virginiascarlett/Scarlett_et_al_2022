#!/bin/bash
#Run this script from the dir ABOVE candidatePseudogeneRegions/
#Need yn00.ctl in the same dir as this script.
#This script is designed to be run in a loop. 
#It will change the first two lines of yn00.ctl to match the current gene,
#then run the yn00 program on that gene.

for fname in $( ls candidatePseudogeneRegions/*.readyForPAML.fasta ); do 
	line="seqfile = $fname"
	#sed -i "1s/.*/$line/" yn00.ctl
	sed -i "1c\\$line" yn00.ctl
	geneID=${fname%.readyForPAML.fasta} 
	line="outfile = $geneID.yn"
	#sed -i "2s/.*/$line/" yn00.ctl
	sed -i "2c\\$line" yn00.ctl
	/opt/paml4.9j/bin/yn00 >> yn00.stdout  #EDIT THIS LINE
done