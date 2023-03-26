#!/usr/bin/python3

#### CODE FOR FINDING UNIQUE-KMER FOR EACH CONTIG IN A MULTI FASTA FILE ######

# module for generating the seq object
from Bio import SeqIO

#######################################################################

# read fasta file from command line argument
fasta_file = "/home/smanavalan/test1.fasta"
# generates a list of seq obj
list_seq_obj = list(SeqIO.parse(fasta_file, 'fasta'))
# path to write the output
output_file = "/home/smanavalan/testresult.fasta"
# kmer size
k_size = 25

#########################################################################

# generate k-mers of user decided size
# number of kmers of size possible from a string of len(seq) is 
# given by the formula len(seq) - k + 1 since the range starts from 0 
# the number of times this loop runs will be equal to number of kmers possible
def generate_kmers(seq, k_size):
    # empty list to store kmers
    kmers = []
    for i in range(len(seq) - k_size + 1):
        # starting postion to the lenght of the kmer
        kmer = seq[i:i + k_size]
        # skip kmers with N
        if 'N' in kmer:
            continue
        # add to kmer list
        kmers.append(kmer)
    return kmers  # return a string of kmer

##################################################################

# generate k-mers for each contig and its rev complement
kmers_dict = {}
for each_seq_obj in list_seq_obj:
    seq = str(each_seq_obj.seq)
    rev_seq = str(each_seq_obj.seq.reverse_complement())
    kmers = generate_kmers(seq, k_size)
    #Store the kmers in a set and add them to the dictionary kmers_dict
    # add them to set type helps in removing duplicates
    kmers_dict[each_seq_obj.id] = set(kmers)
    rev_kmers = generate_kmers(rev_seq, k_size)
    # merge the kmers from rev compliment under the same key 
    # the values in the dict are set any duplicates will be removed
    kmers_dict[each_seq_obj.id].update(rev_kmers)

########################################################################

# find unique k-mers for each contig
# iterate through each seq obj in list , pull out the kmers of each seq obj
# from the dict check if the current kmer is present in dic of other
# seq obj, if not add it to a new dict with the key of current seq obj

# create the final dictionary to store unique kmer
unique_kmers_dict = {}
# iterae through all the seq obj 
for each_seq_obj in list_seq_obj:
    # store kmers of the current seq obj in a set
    unique_kmers = set()
    # Iterating through all the k-mers present in dict for each seq obj
    for current_kmer in kmers_dict[each_seq_obj.id]:
        # set a flag to check if the current kmer is uique  
        unique = True
        # Iterating through all seq obj other than current
        # seq object in input list.
        for other_each_seq_obj in list_seq_obj:
            if other_each_seq_obj.id == each_seq_obj.id:
                continue
            if current_kmer in kmers_dict[other_each_seq_obj.id]:
                unique = False
                # Exit inner loop if kmer present in dict of other seq obj
                break
        if unique:
            # add kmer to set if it is not already present
            unique_kmers.add(current_kmer)
    unique_kmers_dict[each_seq_obj.id] = unique_kmers

################################################################

# write output 
# open the file in write mode
file_obj = open(output_file, 'w')

for each_seq_obj in list_seq_obj:
    # pull out all values associated with the key from dict
    all_kmers_of_id = unique_kmers_dict[each_seq_obj.id]
    # skip the key with no unique-kmer present 
    if len(all_kmers_of_id) == 0:
        continue
    # writes a header line ( contig id) 
    file_obj.write('>' + each_seq_obj.id + '\n')
    # write each kmer present under the key followed by new line
    for each_kmer in all_kmers_of_id:
        file_obj.write(each_kmer + '\n')
file_obj.close()

###################################################################
