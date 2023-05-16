#!/usr/bin/python3

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
# API for MySQL database
import MySQLdb
# for calcualting the melting temperature of sequence
from Bio.SeqUtils import MeltingTemp as mtmp

# argparse library to handle command line argument , the script except
# 5 mandatory arguments > file path , size of kmer   hamming distance and Tm
parser = argparse.ArgumentParser(
    description="Tool for identifying unique oligos from Multi-FASTA record."
)

parser.add_argument(
    "-f",
    required=True,
    help="Path to the multi-FASTA file"
)

parser.add_argument(
    "-k",
    type=int,
    required=True,
    help="Required oligomer size"
)

parser.add_argument(
    "-d",
    type=int,
    required=True,
    help="Hamming distance threshold"
)

parser.add_argument(
    "-min_tm",
    type=float,
    required=True,
    help="Minimum Tm in Celsius"
)

parser.add_argument(
    "-max_tm",
    type=float,
    required=True,
    help="Maximum Tm in Celsius"
)

args = parser.parse_args()

fasta_file = args.f
k_size = args.k
hamming_threshold = args.d
min_tm = args.min_tm
max_tm = args.max_tm

# Establish MySQL database connection using MySQLdb API
# Replace the below paramters according to database configuration.
connection = MySQLdb.connect(
    host="localhost",
    user="smanavalan",
    passwd="",
    db="smanavalan"
)

# Create a cursor object that allows to interact with the database
cursor = connection.cursor()


# func for identifying repeat stretch in the sequence
# returns True if repeat stretch identified
def repeat_stretch(kmer, stretch_length=5):
    for i in range(len(kmer) - stretch_length + 1):
        stretch = kmer[i:i + stretch_length]
        if len(set(stretch)) == 1:
            return True
    return False


# func for calculating Tm based on nearest neighbor thermodynamics.
# returns True if Tm is within the specified range
def tm_range(kmer, min_tm, max_tm):
    tm = mtmp.Tm_NN(kmer)
    return min_tm <= tm <= max_tm


# create an instance of the PairwiseAligner class and configure parameters
# used for checking self annealing properities of kmer using self_dimer()
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5


def self_dimer(kmer, aligner):
    kmer_seq = Seq(kmer)
    # obtain the reverse complement of kmer sequence
    revcomp_seq = kmer_seq.reverse_complement()
    # perform pairwise sequence alignment using Smith-Waterman algorithm
    alignments = aligner.align(kmer_seq, revcomp_seq)
    if alignments:
        best_alignment = alignments[0]
        # calculate the percent identity between the aligned sequences
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
        gc_content = (kmer.count('G') + kmer.count('C')) / k_size * 100
        if gc_content < 40 or gc_content > 60:
            continue
        if not tm_range(kmer, min_tm, max_tm):
            continue
        percent_identity = self_dimer(kmer, aligner)
        # discard kmers with self-complementarity
        if percent_identity > 80:
            continue
        # yield kmer as generator object
        yield kmer


# function to read multi-fasta file and return fasta as generator obj
def read_fasta_file(fasta_file):
    # open mult-fasta file for reading
    with open(fasta_file, "r") as out_file:
        # iterate over the fasta records, output is seq-obj for each fasta
        for fasta in SeqIO.parse(out_file, "fasta"):
            # yield generator obj of type seq-obj
            yield fasta


# generates a list of seq obj for iterting to generate kmers for each obj
list_seq_obj = list(read_fasta_file(fasta_file))

# Create table contig with contig_id as the column, this is the parent
# table storing all the records that are processed, irrespective of
# whether they have unique kmer
cursor.execute(
    "CREATE TABLE IF NOT EXISTS contigs ("
    "contig_id VARCHAR(255) PRIMARY KEY"
    ");"
)

# creates child table unique_kmers with 3 columns > id, contig_id, and kmer
# contig_id > foreign key references > contig_id in contigs table.
# This schema creates a one to many relationship where one contig can have
# many unique kmers. contig_id here will only exist if it has at least one
# unique kmer/oliogs
cursor.execute(
    "CREATE TABLE IF NOT EXISTS unique_kmers ("
    "id INT AUTO_INCREMENT PRIMARY KEY, "
    "contig_id VARCHAR(255), "
    "kmer VARCHAR(255) NOT NULL, "
    "FOREIGN KEY(contig_id) REFERENCES contigs(contig_id)"
    ");"
)

# generate k-mers for each of the record, string generated from both strands
# passed to genearate_kmers().Output is stored in a dict with seq-obj record
# as the key, set() removes any duplicates in the kmers, update() remove
# duplicates between kmers from both the strands for each id
kmers_dict = {}
for each_seq_obj in list_seq_obj:
    seq = str(each_seq_obj.seq)
    rev_seq = str(each_seq_obj.seq.reverse_complement())
    kmers = generate_kmers(seq, k_size)
    kmers_dict[each_seq_obj.id] = set(kmers)
    rev_kmers = generate_kmers(rev_seq, k_size)
    kmers_dict[each_seq_obj.id].update(rev_kmers)

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
            # if not unique stop comparision
            for other_kmer in kmers_dict[other_each_seq_obj.id]:
                if hamming(kmer, other_kmer) <= hamming_threshold:
                    unique = False
                    # if not unique stop iterating through other kmers
                    break

            if not unique:
                # if not unique stop iteration of other-seq-objects and
                # jump to next kmer
                break

        # if unique insert the id, and unique kmer into unique_kmer table
        if unique:
            cursor.execute(
                "INSERT INTO unique_kmers (contig_id, kmer) VALUES (%s, %s);",
                (each_seq_obj.id, kmer)
            )

# commit the insert transction
connection.commit()
# close cursor
cursor.close()
# close the connection to the database
connection.close()
