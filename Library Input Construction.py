import sys
import re

#Open .csv file containing regions upstream of target TSS
#Format of .csv: [pool #, strand, sequence upstream of TSS]
seq_file=open("Human Pos Strand lincRNA UCSC Extractions.txt",'r')
library_file=open("Human Pos Strand lincRNA Promoter Sequences.txt",'w')
strand='+'
pool=0

for extraction in seq_file:
    if extraction[0]==">":
        pool+=1
        m = re.match(r".+=(?P<chrom>\w+):(?P<start>\w+)-(?P<stop>\w+)",extraction)
        elements=[str(pool),m.group('chrom'),m.group('start'),m.group('stop'),strand]
        library_file.write('\n')
        library_file.write(' '.join(elements))
        library_file.write(' ')

    else:
        library_file.write(extraction.rstrip())

#Close files
seq_file.close()
library_file.close()
