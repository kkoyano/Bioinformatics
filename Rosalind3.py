
from collections import Counter
from collections import defaultdict
from operator import itemgetter
from itertools import *
import numpy as np
import sys


KMER_LEN = 4
K = 5


########### OPEN FILE DEFS #################################

def read_reads(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        read = line.strip()

        all_reads.append(read)
    return all_reads

def openfile(fn):
    f = open(fn, 'r')
    f.readline() #skip first header
    string = ''
    for line in f:
        string+=line.strip()

    return string

def openfile_edit_distance(fn):
    f = open(fn, 'r')
    f.readline() #skip first header
    string1 = ''
    string2= ''
    for line in f:
        if line.startswith('>'):
            break
        string1+=line.strip()

    for line in f:
        string2+=line.strip()

def open_superstring(fn):
    f = open(fn, 'r')
    f.readline() #skip first header

    dna_strings = []
    sequence = ''
    for line in f:

        if line.startswith('>'):
            dna_strings.append(sequence)
            sequence = ''
            continue
        sequence+=(line.strip())
    dna_strings.append(sequence) #capture the last sequence
    # print len(dna_strings)
    # print dna_strings
    return dna_strings
#######################################################################

def simple_de_bruijn(sequence_reads):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of
size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn
graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the
set

    NOTE: THIS CONTAINS THE INFO FOR ROSALIND PROBLEMS 10 AND 11

    """

    de_bruijn_graph = defaultdict(list)
    output_fn = open('pattern.out.txt', 'w')
    for read in sequence_reads:

        k= len(read) -1
        key = read[:k]
        value= read[1:]
      
      
      
        # de_bruijn_graph[key]= value
        de_bruijn_graph[key].append(value)
        # print de_bruijn_graph
    for key in sorted(de_bruijn_graph):
        key_out= '{} -> '.format(key)

        values = ''
        count = 0
        for value in de_bruijn_graph[key]:
            if count >= 1 and count < len(value):
                values+= ', '
            values+= value
            count +=1

        out = key_out + values
        # print out
        output_fn.write(out)
        output_fn.write('\n')

        #unblock this for problem 10 and also add default counter instead of just a regular list. basically use michaels code
        # kmers = [read[i: i + k] for i in range(len(read) - k)]
        # for i in range(len(kmers) - 1):
        #     pvs_kmer = kmers[i]
        #     next_kmer = kmers[i + 1]
        #     de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen
    # enough.
    # de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if
# de_bruijn_counter[key][val] > 0}
#                        for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
#     de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if
# de_bruijn_graph[key]}
#     # print de_bruijn_graph
#     print len(de_bruijn_graph)
    return de_bruijn_graph


def de_bruijn_reassemble(de_bruijn_graph):
    """
    Traverses the DeBruijn Graph created by simple_de_bruijn and
    returns contigs that come from it.
    :param de_bruijn_graph: A De Bruijn Graph
    :return: a list of the
    """
    # print de_bruijn_graph
    assembled_strings = []
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        if n_values == 0:
            break
        good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]
        # print good_starts
        # You may want to find a better start
        # position by looking at in and out-degrees,
        # but this will work okay.
        current_point = good_starts[0]
        assembled_string = current_point
        # current_point = good_starts[0]
        # assembled_string = current_point
        while True:
            try:
                next_values = de_bruijn_graph[current_point] #edge from current point
                next_edge = next_values.pop()
                assembled_string += next_edge[-1]
                de_bruijn_graph[current_point] = next_values
                current_point = next_edge
            except KeyError:
                assembled_strings.append(assembled_string)
                break
    return assembled_strings


def db_of_string(fn):
    f = open(fn, 'r')
    k = f.readline()
    k = k.strip()
    k = int(k)
    string = f.readline()
    string = string.strip()

    # kmer_list = [string[i:i+k] for i in range(0,(len(string)-k))]
    # return kmer_list
    de_bruijn_graph = defaultdict(list)
    k= k-1
    output_fn = open('db_from_string.out.txt', 'w')
    for i in range(0,(len(string)-k)):
        # print read

        key = string[i:i+k]
        value= string[i+1:i+k+1]

        # de_bruijn_graph[key]= value
        de_bruijn_graph[key].append(value)
        # print de_bruijn_graph
    for key in sorted(de_bruijn_graph):
        key_out= '{} -> '.format(key)

        values = ''
        count = 0
        for value in de_bruijn_graph[key]:
            if count >= 1 and count < len(value):
                values+= ', '
            values+= value
            count +=1

        out = key_out + values
        # print out
        output_fn.write(out)
        output_fn.write('\n')

if __name__ == "__main__":
    # input= sys.argv[1]
    fn = sys.argv[1]


    db = db_of_string(fn)
    # print kmer_list
    # reads = read_reads(fn)
    # db_graph = simple_de_bruijn(kmer_list)
    # db_graph= db(reads)
    # print db_graph
    # out = de_bruijn_reassemble(db_graph)




