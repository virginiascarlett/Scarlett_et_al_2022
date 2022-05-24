Here's how I ran TEMP2 to get B. hybridum TE polymorphisms.
First, I got the TE library by running my TE pipeline. Then I
annotated the reference genome with this library using RepeatMasker. 

To get TE polymorphisms among natural B. hybridum lines, I used ABR113 as
my reference genome. 

In my dir /global/cscratch1/sd/vstartag/jacob_internship/McClintock_inputs,
I ran prep_McClintock_raw_gff.py. That script concatenates the gffs 
RepeatMasker produced for each chromosome (since I ran RepeatMasker one 
chromosome at a time). It also depends on those being in a particular file 
structure described in the config file for my TE pipeline, timeForTE.py.
That produces a TE gff, e.g. ABR113_RepeatMasker_TElib.gff

Then I ran prep_McClintock_final_inputs.py. This adjusts the TE library 
headers to be compatible with McClintock/TEMP2, which literally doesn't
accept any special characters at all (very annoying).
Need to specify the genome and then the library filename like so e.g.
python3 prep_McClintock_final_inputs.py ABR113 TElib.ABR113.fa 
This script will create a new TE library, e.g. TElib.ABR113.clean.fa, 
and a new TE gff e.g. ABR113_RepeatMasker_TElib.clean.gff.
This script will also produce a TSV file with family classifications, e.g.
ABR113_TElib_McClintockfams.tsv. Then I ran mcclintock with a job array as in
/global/cscratch1/sd/vstartag/jacob_internship/Bhybridum/mcclint.sh

To run on cori, use this conda environment:
module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/mcclintock

