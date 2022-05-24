
FILES=*.trimmed
for f in $FILES
do
	IFS='.' 
	read -ra myarray <<<"$f"
	basename="${myarray[0]}"
	distmat "${basename}.trimmed" -nucmethod 2 -outfile "${basename}.distmat" 
done