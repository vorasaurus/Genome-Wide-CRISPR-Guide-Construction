import sys
import re

#Open .csv file containing regions upstream of target TSS
#Format of .csv: [pool #, strand, sequence upstream of TSS]
seq_file=open("Figure 2 Library Targets 6-6-14.csv",'r')
target_file=open("Activator Library Targets 6-6-14.txt",'w')

#Generate translation table to generate compliment sequence
compliment = str.maketrans('ATCG','TAGC')

#Function to search for matches to pam on strand, starting from TSS
def search_pam(pam_pattern,seq):

    #Search for start positions of pam sites in sequence, starting from TSS
    starts = [match.start() for match in re.finditer(pam_pattern, seq)]

    #Return start positions for target sites identified by pam pattern recognition
    return starts

#Function to select target sequences given target strand, search strand, and start positions
def generate_targets(pool,strand,seq,starts_pos,starts_neg):

    #At most, select first 15 targets in list of positions
    if len(starts_pos)>= 15:
        starts_pos = starts_pos[0:15] 
    if len(starts_neg)>=15:
        starts_neg = starts_neg[0:15]

    #Select target sequences, given target is on + strand
    if strand == '+':

        #For pam start positions on positive strand
        for position in starts_pos:
            #select 20 bp upstream of pam and de-reverse seq. to get original seq.
            target=seq[position+3:position+23][::-1]
            #Generate target output to the target file
            position=''.join(['-',str(position)])
            search_strand='+'
            generate_target_output(pool,target,search_strand,position)

        #For pam start positions on negative strand
        for position in starts_neg:
            #select 20 bp upstream of pam and compliment seq. to get reverse comp.
            target=seq[position:position+20]
            target=target.translate(compliment)
            #Generate target output to target file
            position=''.join(['-',str(position)])
            search_strand='-'
            generate_target_output(pool,target,search_strand,position)

    #Select target sequences, given target is on - strand
    elif strand == '-':
        #For pam start positions on positive strand
        for position in starts_pos:
            #select 20 bp upstream of pam 
            target=seq[position:position+20]
            #Generate target output to the target file
            position=''.join(['+',str(position)])
            search_strand='+'
            generate_target_output(pool,target,search_strand,position)

        #For pam start positions on negative strand
        for position in starts_neg:
            #select 20 bp upstream of pam and reverse comp.
            target=seq[position+3:position+23][::-1]
            target=target.translate(compliment)
            #Generate target output to target file
            position=''.join(['+',str(position)])
            search_strand='-'
            generate_target_output(pool,target,search_strand,position)

#Function to write target information [pool#, target seq., strand, position from TSS]
def generate_target_output(pool,target,strand,position):
    #Compile target information into list
    target_info=[str(pool),target,strand,position]
    #Splice together with commas
    target_line=','.join(target_info)
    #Write to target file
    #target_file.write(target_line)
    target_file.write(target_line)
    target_file.write('\n')

#Start from TSS. Generate 15 spacer sequences, separated by >13bp, both strands.

#Read lines of locus information
for locus_info in seq_file:
    
    #Parse locus information to get pool#, + or - strand, and TSS upstream seq
    pool,target_strand,seq=locus_info.split(',')

    #Set regex search pattern for pam according to target and search strand
    if target_strand=='+':
        #Search backwards for targets on + strand, in order to start from TSS
        seq = seq[::-1]
        #Search pos strand for GGNN{20}N{13} with 13 bp spacing between guides
        pam_pos_strand = 'GG.{34}'
        #Search neg strand for N{20}NCCN{13} with 13 bp spacing between guides
        pam_neg_strand = '.{21}CC.{13}'

    elif target_strand == '-':
        #Search pos strand for N{20}NGGN{13} with 13 bp spacing between guides
        pam_pos_strand = '.{21}GG.{13}'
        #Search neg strand for CCNN{20}N{13} with 13 bp spacing between guides
        pam_neg_strand = 'CC.{34}'

    #Calculate start positions in positive strand for instances of pam recognition sequence
    starts_pos_strand = search_pam(pam_pos_strand,seq)
    
    #Calculate start positions in positive strand for instances of pam recognition
    starts_neg_strand = search_pam(pam_neg_strand,seq)
    
    #Generate target sequences given start positions on both + and - strands
    generate_targets(pool,target_strand,seq,starts_pos_strand,starts_neg_strand)

#Close files
seq_file.close()
target_file.close()
