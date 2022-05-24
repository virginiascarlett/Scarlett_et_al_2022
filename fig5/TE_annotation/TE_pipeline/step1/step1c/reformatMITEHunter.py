"""
A script to change ugly MITE-Hunter headers (contain spaces, special symbols, etc...)
into simple headers ready for RepeatMasker. No args needed.
"""
import timeForTE

MITEfile = open('MITE-Hunter.allresults.fa', 'r')

counter = 0
lines = []
for line in MITEfile.readlines():
    if line.startswith('>'):
        counter += 1
        lines.append('>mite{}#DNA/Unknown\n'.format(str(counter)))
    else:
        lines.append(line)
MITEfile.close()
        
with open('tempfile', 'w') as outF:
    for line in lines:
        outF.write(line)
        
timeForTE.reformat_fasta('tempfile', 'MITE-Hunter.allresults.fa.RMready', 50)
