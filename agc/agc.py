#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    Reads a compressed fasta file an returns sequnces long enough

    Parameters
    ----------
    amplicon_file: str
    Path to the file
    minseqlen: int
    the minimum length of a read sequence

    Returns
    -------
    A generator containing all sequences passing the cut
    """
    # Only works for single-line fastas
    with gzip.open(amplicon_file, "rt") as filein:
        line = filein.readline().rstrip()
        seq = ""
        while line:
            if line.startswith(">"):
                if len(seq) > minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line
            line = filein.readline().rstrip()
        if len(seq) > minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Creates a list of lists of sequence, count from a file containing sequences

    Parameters
    ----------
    amplicon_file: str
    pat to the file
    minseqlen: int
    minimym length of sequences
    mincount: int
    minimum count of sequences

    Returns
    -------
    A generator of the [sequence, count] lists
    """
    dic_count = {seq:count for seq, count in Counter(read_fasta(amplicon_file, minseqlen)).items() if count >= mincount}
    for seq, count in sorted(dic_count.items(), key=lambda x:x[1], reverse=True):
        yield [seq, count]



def get_identity(alignment_list):
    """Returns the identity scores from aligned sequences in a list

    Parameters
    ----------
    alignment_list: list
    list of lists of aligned sequences

    Returns
    -------
    A list of alignement scores
    """
    return sum([1 for i in range(len(alignment_list[0])) if alignment_list[0][i] == alignment_list[1][i]])/len(alignment_list[0]) * 100


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size=None, kmer_size=None):
    """
    Calculates an OTU clustering based on a fasta file.

    Parameters
    ----------
    amplicon_file: str
    pat to the file
    minseqlen: int
    minimym length of sequences
    mincount: int
    minimum count of sequences

    Returns
    -------
    A list of lists of OTU sequences and their respective counts
    """
    lists_of_derep = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    dict_otus = dict([lists_of_derep[0]])
    for seq, count in lists_of_derep[1:]:
        refseq = ""
        for seq_o in dict_otus:
            align = nw.global_align(
                seq_o,
                seq,
                gap_open=-1,
                gap_extend=-1,
                matrix=os.path.abspath(os.path.join(os.path.dirname(__name__),"agc/MATCH"))
            )
            if get_identity(align) >= 97:
                refseq = seq_o
        if refseq:
            dict_otus[refseq] += count
        else:
            dict_otus[seq] = count
    return list(dict_otus.items())





def write_OTU(OTU_list, output_file):
    """
    Writes the given OTU_list to an output file

    Parameters
    ----------
    OTU_list: list
    List of tuples of OTUs sequences and their counts
    output_file: str
    Path to an output file
    """
    with open(output_file, "w") as fileout:
        for i in range(len(OTU_list)):
            fileout.write(f">OTU_{i+1} occurrence:{OTU_list[i][1]}\n")
            fileout.write(textwrap.fill(OTU_list[i][0], 80) + "\n")


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici

#==============================================================
# Chimera removal section
#==============================================================

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def get_chunks(sequence, chunk_size):
    """Split sequences in a least 4 chunks
    """
    pass

def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass

def detect_chimera(perc_identity_matrix):
    pass

def search_mates(kmer_dict, sequence, kmer_size):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


if __name__ == '__main__':
    main()
