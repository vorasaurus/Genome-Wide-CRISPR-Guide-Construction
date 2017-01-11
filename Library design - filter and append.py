import sys
import re

#Open .txt file containing library target sites
#Format of .csv: [pool #, target sequence, strand, distance from TSS]
target_file=open("Mouse lincRNA Library 1.3.txt",'r')
filter_file=open("Mouse lincRNA Library 1.3 filtered and appended.txt",'w')


for target in target_file:
    pool, protospacer, strand, distance = target.split(',')
    N_filter = re.search('N', protospacer)
    BsmBI_filter_one = re.search('CGTCTC', protospacer)
    BsmBI_filter_two = re.search('GAGACG', protospacer)
    if not N_filter and not BsmBI_filter_one and not BsmBI_filter_two:
        new_protospacer = 'CGTCTCACACCG'+protospacer+'GTTTTGAGACG'
        filtered_info=[str(pool),new_protospacer,strand,distance]
        filtered_line=','.join(filtered_info)
        filter_file.write(filtered_line)


#Close files
target_file.close()
filter_file.close()
