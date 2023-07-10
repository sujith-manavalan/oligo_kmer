
"""
Code for identification of unique kmer for printing in a micro-array chip
STEP 1 >>>>>>>>>>>
# reads multi-FASTA file using SeqIO.parse()
# Each sequence in the record, along with reverse complement is processed to
generate all possible kmers of a given size.
STEP 2 >>>>>>>>>>>>
Each kmer undergo filtering:
    # should not contain'N' bases.
    # should not contain homopolymer stretches of length 5 or more.
    # should be within a specified GC range.
    # Tm should be within a specified range.
    # should not self-dimerize( checked using reverse smith-watermna algorithm)
STEP 3 >>>>>>>>>>>>>
 # Unique kmers identified by comparing each kmer of a record against all
 kmers from all other records(slowest part of the code).
 # kmer is unique if no match with any kmers (either perfect string matching
 or partial matching(hamming))
STEP 4 >>>>>>>>>>>>>>>
unique kmers written to multi-FASTA file or stored in a MySQL database.
MySQL > creates two tables, contigs and unique_kmers
contigs table stores the IDs of all processed sequnece, and unique_kmers table
stores the unique kmers identified, along with the ID of the sequence(one to
many schema)

"""

# for handling command line argument
import argparse
# I/O interface for biopython
from Bio import SeqIO
# used for local alignments between two sequences
from Bio.Align import PairwiseAligner
# for  creating Seq class objects
from Bio.Seq import Seq
# for calculating Hamming distance between two sequences
from distance import hamming
# for calcualting the melting temperature of sequence
from Bio.SeqUtils import MeltingTemp as mtmp
# for multi-processing (built in library)
# imports the Pool class from the multiprocessing module
from multiprocessing import Pool
# API for MySQL database
import MySQLdb

# argparse library to handle command line argument , the script except
# 5 mandatory arguments > file path , size of kmer   hamming distance and Tm
# parser = argparse.ArgumentParser(
#     description="Tool for identifying unique oligos from Multi-FASTA record.",
#     formatter_class=argparse.RawTextHelpFormatter
# )

# parser.add_argument(
#     "-f",
#     required=True,
#     help="Path to the multi-FASTA file"
# )

# parser.add_argument(
#     "-k",
#     type=int,
#     required=True,
#     help="Required oligomer size"
# )

# parser.add_argument(
#     "-d",
#     type=int,
#     help="Hamming distance threshold"
# )

# parser.add_argument(
#     "-min_tm",
#     type=float,
#     help="Minimum Tm in Celsius"
# )

# parser.add_argument(
#     "-max_tm",
#     type=float,
#     help="Maximum Tm in Celsius"
# )
# parser.add_argument(
#     "-min_gc",
#     type=float,
#     help="Minimum GC content"
# )

# parser.add_argument(
#     "-max_gc",
#     type=float,
#     help="Maximum GC content"
# )
# parser.add_argument(
#     "-Na",
#     type=float,
#     help="Molar Na concentration(for salt adjusted Tm calculation)"
# )
# parser.add_argument(
#     "-i",
#     type=int,
#     help="Sequence similarity percentage threshold(Using Smith-Waterman\n"
#     "algorithm to remove oligos with self-dimerization tendency)"
# )
# parser.add_argument(
#     "-o",
#     help="Optional path to the output multi-FASTA file"
# )

# args = parser.parse_args()

fasta_file = "mini.fasta"  # input file
hamming_threshold = None   # pass None if hamming is not required
k_size = 24
min_tm = 1
max_tm = 99
min_gc = 40
max_gc = 60
Na = None
input_identity = 80
output_file = "test_output.fa"


# use MySQL database connection using MySQLdb API only if not output file is
# mentioned, replace the below paramters according to database configuration.
if not output_file:
    connection = MySQLdb.connect(
        host="localhost",
        user="xxxx",
        passwd="xxxx",
        db="suj"
    )

    # Create a cursor object that allows to interact with the database
    cursor = connection.cursor()

    # Create table named contigs with contig_id as the column, this is the
    # parent table storing all the records that are processed, irrespective of
    # whether they have unique kmer
    cursor.execute(
        "CREATE TABLE IF NOT EXISTS contigs ("
        "contig_id VARCHAR(255) PRIMARY KEY"
        ");"
    )

    # creates child table unique_kmers with 3 columns > id, contig_id, and kmer
    # contig_id > foreign key references > contig_id in contigs table.
    # schema creates a one to many relationship where one contig_id can
    # have many unique kmers. contig_id here will only exist if it has at
    # least one unique kmer/oliogs. Schema also helps in retreival of records
    # with no unique kmer
    cursor.execute(
        "CREATE TABLE IF NOT EXISTS unique_kmers ("
        "id INT AUTO_INCREMENT PRIMARY KEY, "
        "contig_id VARCHAR(255), "
        "kmer VARCHAR(255) NOT NULL, "
        "FOREIGN KEY(contig_id) REFERENCES contigs(contig_id)"
        ");"
    )


# func for identifying repeat stretch of a nucleotide(of length 5 or more) in
# the sequence, returns True if repeat stretch(homopolymer) identified
def repeat_stretch(kmer, stretch_length=5):
    # loop for iterating thorugh substring of length = stretch_length
    for i in range(len(kmer) - stretch_length + 1):
        # get the substring from kmer
        stretch = kmer[i:i + stretch_length]
        # using set() check if all nucleotides are same
        if len(set(stretch)) == 1:
            return True
    return False


# func for calculating Tm based on nearest neighbor thermodynamics.
# returns True if Tm is within the specified range. If Na conc is given
# the  calculates Tm using salt adjusted formula
def tm_range(kmer, min_tm, max_tm):
    if Na is not None:
        tm = mtmp.Tm_NN(kmer, Na=Na)
    else:
        tm = mtmp.Tm_NN(kmer)
    return min_tm <= tm <= max_tm


# create an object using PairwiseAligner class and configure attributes
# used for checking self annealing properities of kmer using self_dimer()
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5


# func for checking whether kmer undergo self-dimerization, Smith
# water algorithm(parameters above) is used here for checking local
# alignment score. func takes kmer and aligner object
def self_dimer(kmer, aligner):
    # covert to seq object to use reverse_complement()
    kmer_seq = Seq(kmer)
    # obtain the reverse complement of kmer sequence
    revcomp_seq = kmer_seq.reverse_complement()
    # perform pairwise sequence alignment using Smith-Waterman algorithm
    alignments = aligner.align(kmer_seq, revcomp_seq)
    # if alignment exsists retreive the best alignment[0]
    if alignments:
        best_alignment = alignments[0]
        # Calculate the alignment score as a percentage of
        # maximum possible score
        percent_identity = (best_alignment.score / len(kmer)) * 100
    # if no aligment condition then set the value to 0
    else:
        percent_identity = 0

    return percent_identity


# function to generate kmers of user decided size, number of kmers of k_size
# possible from a string of len(seq) is given by the formula
# len(seq) - k_size + 1 since the range starts from 0
# gc_content, repeate_stretch, self-complimentary and kmers not in Tm range
# are fitered out in the function
def generate_kmers(seq, k_size):
    for i in range(len(seq) - k_size + 1):
        kmer = seq[i:i + k_size]
        if 'N' in kmer:
            continue
        if repeat_stretch(kmer):
            continue
        if min_gc is not None and max_gc is not None:
            gc_content = (kmer.count('G') + kmer.count('C')) / k_size * 100
            if gc_content < min_gc or gc_content > max_gc:
                continue
        # func call only when user has supplied the input
        if min_tm is not None and max_tm is not None:
            if not tm_range(kmer, min_tm, max_tm):
                continue
        if input_identity is not None:
            percent_identity = self_dimer(kmer, aligner)
            # discard kmers with self-complementarity
            if percent_identity > 80:
                continue
        # yield kmer as generator object
        yield kmer


# function to read multi-fasta file and return fasta as generator obj
def read_fasta_file(fasta_file):
    # open mult-fasta file for reading with out_file as the file object
    with open(fasta_file, "r") as file_obj:
        # iterate over the fasta records, output is seq-obj for each fasta
        for fasta in SeqIO.parse(file_obj, "fasta"):
            # yield generator obj of type seq-obj
            yield fasta


# generates a list of seq obj for itertion to generate kmers as well as to
# check uniquness ( avoid multiple reading)
list_seq_obj = list(read_fasta_file(fasta_file))


# function to generate filtered kmers from both the strands, requires a seq obj
# and the kmer size. Genrate string from seq obj then calls generate_kmers()
# and use set() to remove duplicates in kmers . returns a tuple containing
# (seq_obj.id, kmer1 , kmer2). the returned tuple is iterated to
# generate kmer_dict{}. This function will be prallellized
def parallel_job(seq_obj, k_size):
    seq = str(seq_obj.seq)
    rev_seq = str(seq_obj.seq.reverse_complement())
    return (seq_obj.id, set(generate_kmers(seq, k_size)),
            set(generate_kmers(rev_seq, k_size)))


# create an instance of multiprocessing pool class, all cpus will be used
pool = Pool()

# generate arguments for the parallel_job(), this is a list of tuples[(),()]
# where each tuple contains a seq obj and the kmer size
parallel_job_arguments = []
for seq_obj in list_seq_obj:
    parallel_job_arguments.append((seq_obj, k_size))
# Use starmap() of multiprocessing Pool to apply the parallel_job
# function to each set of arguments. starmap() unpacks each tuple in the list
#  and passes the elements as separate arguments to parallel_job()
all_kmers = pool.starmap(parallel_job, parallel_job_arguments)

# close the multiprocessing pool
pool.close()

print("Finished generating kmers.", flush=True)

# dict to store kmers
kmers_dict = {}
# Each result is tuple containing ID and 2 sets of kmers(compli and rev)
# results will still be in the order of the input iterable that is
# parallel_job_arguments
for each_result in all_kmers:
    id = each_result[0]
    kmers = each_result[1]
    rev_kmers = each_result[2]
    kmers_dict[id] = kmers.union(rev_kmers)

############################################################################

if output_file:
    # this block is for writing the output to multi-fasta file
    with open(output_file, "w") as out_fasta:
        # find unique k-mers for each contig by comparing kmers of each seq-obj
        # against kmers of all other seq-obj using hamming distance
        for each_seq_obj in list_seq_obj:
            # Insert id of seq-obj into fasta
            out_fasta.write('>' + each_seq_obj.id + '\n')

            # Iterating through all k-mers present in dict for each seq-obj
            for kmer in kmers_dict[each_seq_obj.id]:
                # flag true as long as the kmer is unique
                unique = True

                # Iterating over all seq objects
                for other_each_seq_obj in list_seq_obj:
                    # avoid comparsion of same seq-object
                    if other_each_seq_obj.id == each_seq_obj.id:
                        continue

                    # iterating over all kmers of other seq-obj and comparing
                    # current kmer with other kmers using hamming distance,
                    # if not unique stop comparision, according to user input
                    # uses either hamming or perfect string matching
                    for other_kmer in kmers_dict[other_each_seq_obj.id]:
                        if hamming_threshold is None:
                            if kmer == other_kmer:
                                unique = False
                                # if not unique stop iterating through
                                # other kmers
                                break
                        else:
                            if hamming(kmer, other_kmer) <= hamming_threshold:
                                unique = False
                                # if not unique stop iterating through
                                # other kmers
                                break

                    if not unique:
                        # if not unique stop iteration of other-seq-objects and
                        # jump to next kmer
                        break

                # if unique insert the id, and unique kmer into
                # unique_kmer table
                if unique:
                    out_fasta.write(kmer + '\n')

else:
    # this block is for SQL databse
    # find unique k-mers for each contig by comparing kmers of each seq-obj
    # against kmers of all other seq-obj using hamming distance
    for each_seq_obj in list_seq_obj:
        # Insert id of seq-obj into MySQL parent contigs table
        cursor.execute(
            "INSERT INTO contigs (contig_id) VALUES (%s);", (each_seq_obj.id,)
                    )

        # Iterating through all k-mers present in dict for each seq-obj
        for kmer in kmers_dict[each_seq_obj.id]:
            # flag true as long as the kmer is unique
            unique = True

            # Iterating over all seq objects
            for other_each_seq_obj in list_seq_obj:
                # avoid comparsion of same seq-object
                if other_each_seq_obj.id == each_seq_obj.id:
                    continue

                # iterating over all kmers of other seq-obj and comparing
                # current kmer with other kmers using hamming distance,
                # if not unique stop comparision, depending on the user input
                # uses either hamming or perfect string matching
                for other_kmer in kmers_dict[other_each_seq_obj.id]:
                    if hamming_threshold is None:
                        if kmer == other_kmer:
                            unique = False
                            # if not unique stop iterating through other kmers
                            break
                    else:
                        if hamming(kmer, other_kmer) <= hamming_threshold:
                            unique = False
                            # if not unique stop iterating through other kmers
                            break

                if not unique:
                    # if not unique stop iteration of other-seq-objects and
                    # jump to next kmer
                    break

            # if unique insert the id, and unique kmer into unique_kmer table
            # in the column named contig_id and kmer
            if unique:
                cursor.execute(
                    "INSERT INTO unique_kmers (contig_id, kmer) "
                    "VALUES (%s, %s);",
                    (each_seq_obj.id, kmer)
                )

    # commit the insert transction
    connection.commit()
    # close cursor
    cursor.close()
    # close the connection to the database
    connection.close()

print("Finished writing output file.")
