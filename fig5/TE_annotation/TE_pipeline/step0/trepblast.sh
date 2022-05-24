cd TREP_exemplars_to_blast

FILES=*_first.fa
for f in $FILES
do
	IFS='_' 
	read -ra myarray <<<"$f"
	basename="${myarray[0]}"
	blastn -query "${basename}_first.fa" -subject "${basename}_second.fa" -out "${basename}.aln"
done