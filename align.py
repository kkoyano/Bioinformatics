#!/usr/bin/python
import sys
import numpy as np
from os.path import join
import time



#function that opens a file called read_fn, and puts each read
#(which is separated by each line, or a paired end read is separated by a comma)
# into a list. Each index of the list contains a single read.
def read_reads(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip() #removes white space and new line character
        paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
        #eg: paired ends will be ['AGCTC', 'GAGCT']
        all_reads.append(paired_end_reads)
    return all_reads

#function that takes the reference file, called ref_fn, and concatenates all of the
#reads into a single string, which will be returned as output_reference
def read_reference(ref_fn):
    f = open(ref_fn, 'r')
    first_line = True
    output_reference = ''
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        output_reference += line  # We append each line to the output reference string.
    return output_reference



#########################################    HASHING        ###################################

#################################################################################

def PartitionReference(reference_seq, k=16):
    """
    :param reference_seq: A string which is the reference genome
    :param k: The length of each segment *i.e. break reference into 16nt intervals
    :return: A dictionary of the reference sequence partitioned into keys that are kmers of k length. The value for
    each key is the index at which that kmer first appears.
    
    NOTE: 16 is found from...
    number of max_mismatches allowed = 2 
    split sequence into 3 parts (aka number of mismatches allowed + 1, such that there will be one part that perfectly aligns if it 
    is a valid sequence alignment. 
    read length = 50 bases
    
    k= read length/ (mismatches +1) 
    k= 50/(2+1)
    k= 16.6 
    round to 16. 
    
    
    """
    reference_dict = {}
    i=0
    #initialize each kmer to be assigned to a list--> this is for the case of multiple indexes
    for i in range(len(reference_seq)-k):
        kmer = reference_seq[i:i + k]
        reference_dict[kmer] = []

    #append the index for each time the kmer appears
    for i in range(len(reference_seq)-k):
        kmer = reference_seq[i:i + k]
        reference_dict[kmer].append(i)

    return reference_dict

def PartitionRead(read, k=16):
    """
    :param read: A string which is the read you want to segment
    :param k: the length of nucleotides you want the read to be partitioned (should be same as the k variable in PartitionReference)
    :return: splitsequences, a list that has the read partitioned 
    """

    splitsequences = []
    numsegments = int(round(len(read)/k))
    for i in range(1, numsegments + 1):
        if i==1:
         splitsequences.append(read[0:k])
        else:
            splitsequences.append(read[k*(i-1):k*i])

    return splitsequences #ex: ['CTAC', 'GATC', 'GACT']


def trivial_algorithm(paired_end_reads, ref):
    """

    :param paired_end_reads: Paired-end reads generated from read_reads
    :param ref: A reference genome generated from read_reference
    :return: 2 lists:
                1) a list of alignment locations for each read (all_alignment_locations).
                    The list gives the starting position of the minimum-mismatch alignment of both reads.
                2) a list of the paired-end reads set so that both reads are in their optimal orientation
                   with respect to the reference genome.
    """
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    start = time.clock()

    max_mismatch=2
    ref_dict = PartitionReference(ref) #default segment_length = 25nt
#    print ref_dict
#    print 'GENOME = %d' %len(ref)
#    ref_dict= PartitionReference(ref, segment_length)
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 20 == 0:
            time_passed = (time.clock()-start)/60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print 'Approximately {:.2} minutes remaining'.format(remaining_time)

        for read in read_pair:

            #initialize match score and location
            top_match_score = -1
            top_match_location = -1
            top_match_read = -1
            segment_length = int(round(len(read) /(max_mismatch + 1))) #5
            read_segments= PartitionRead(read,segment_length)
#            print read_segments
            for i in range(max_mismatch+1):
        #        print "this is i = %d" %i
                start = i*segment_length #i=0, start =0, end = 5
                end= min((i+1)*segment_length, len(read))


                matching= [ref_dict[f] for f in read_segments if f in ref_dict]
                if i >= len(matching):
                    break
                segment_match_indexes = matching[i]

                #need a list of locations where the segments match the reference
                for index in segment_match_indexes: #ex: [0, 23, 190]
#                    print 'index %d' %index
                    if index < start or index-start+len(read) > len(ref): #makes sure that index is in bounds
                        continue
                    mismatches = 0
                    score = 0
                    for j in range(0,start): #compare nucleotides before segment
                        score += 1
                        if not read[j] ==ref[index-start+j]:
                            mismatches += 1
                        if mismatches > max_mismatch:
                                break
                    for j in range(end,len(read)): #compare nucleotides after segment
                        score +=1
                        if not read[j] ==ref[index-start+j]:
                            mismatches += 1
                        if mismatches > max_mismatch:
                                break
#                    print 'score = %d' %score
        #maybe createa tuple with the first index being the score and the next being the index of the reference which that score occured,
                    if mismatches <= max_mismatch:
                        if score > top_match_score:
                            top_match_score = score
                            top_match_read = read
                            top_match_location = index-start

######################         If Reverse Strand               ##################
            reversed_read= read[::-1]
            read_segments= PartitionRead(reversed_read,segment_length)
            for i in range(max_mismatch+1):
                matching= [ref_dict[f] for f in read_segments if f in ref_dict]
                if i >= len(matching):
                    break
                segment_match_indexes = matching[i]

                #need a list of locations where the segments match the reference
                for index in segment_match_indexes: #ex: [0, 5, 10]
#                    print 'index %d' %index
                    if index < start or index-start+len(read) > len(ref): #makes sure that index is in bounds
                        continue
                    mismatches = 0
                    score = 0
                    for j in range(0,start): #compare nucleotides before segment
                        score += 1
                        if not reversed_read[j] ==ref[index-start+j]:
                            mismatches += 1
                            if mismatches > max_mismatch:
                                break
                    for j in range(end,len(read)): #compare nucleotides after segment
                        score +=1
                        if not reversed_read[j] ==ref[index-start+j]:
                            mismatches += 1
                            if mismatches > max_mismatch:
                                break
                    if mismatches <= max_mismatch:
                        if score > top_match_score:
                            top_match_score = score
                            top_match_read= reversed_read
                            top_match_location = index-start


#                if top_match_score == -1:
#                    read = ''


            read_alignment_locations.append(top_match_location) #position that obtained the highest score is stored in a list
            output_read_pair.append(top_match_read) #corresponding read is also stored in this list, so indexes match

        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)

    return all_read_alignment_locations, output_read_pairs


def pretty_print_aligned_reads_with_ref(genome_oriented_reads, read_alignments, ref):
    """
    :param genome_oriented_reads: oriented reads generated by trivial_algorithm
    :param read_alignments: alignments generated from trivial_algorithm
    :param ref: reference generated by read_ref
    :return: Returns nothing, but prints the reads aligned to the genome to
     show you what pileup actually *LOOKS* like. You should be able to call SNPs
     by eyeballing the output. However, there are some reads that will not align.
     In the future you'll want to re-check why these reads aren't aligning--the cause
     is usually a structural variation, like an insertion or deletion.
    """

    output_str = ''
    good_alignments = [120 < x[1] - x[0] < 180 for x in read_alignments]
    # There should be 50 + x (90 < x < 110) p between the reads, and we give a little
    # extra space in case there's been a deletion or insertion.  Depending on the type of
    # deletions/insertions
    print genome_oriented_reads[-1]
    print read_alignments[-1]
    best_reads = [genome_oriented_reads[i] for i in range(len(good_alignments))
                  if good_alignments[i]]
    # Remove the reads that do not have a good alignment, or a good reverse alignment.
    best_alignments = [read_alignments[i] for i in range(len(read_alignments))
                       if good_alignments[i]]
    # Take their corresponding alignments
    print len(best_reads)
    print best_alignments[-1]

    aligned_reads = [str(best_reads[i][0]) + '.' * (best_alignments[i][1] - best_alignments[i][0] - 50)
                     + str(best_reads[i][1]) for i in range(len(best_reads))]


    # This turns the reads into strings oriented towards the genome.
    # We get the first read, followed by the correct number of dots to join the first and second reads,
    # and then the second read.

    first_alignment = [x[0] for x in best_alignments]
    alignment_indices = np.argsort(first_alignment)
    sorted_reads = [aligned_reads[i] for i in alignment_indices]
    sorted_alignments = [best_alignments[i] for i in alignment_indices]


    # You don't need to worry too much about how the code block below works--its job is to make it so
    # that a read that starts printing in the third row will continue printing in the third row of the
    # next set of lines.
    active_reads = []
    line_length = 100
    output_str += '\n\n' + '-' * (line_length + 6) + '\n\n'
    for i in range(len(ref) / line_length):
        next_ref = ref[i * line_length: (i + 1) * line_length]
        new_read_indices = [j for j in range(len(sorted_reads))
                            if i * line_length <= sorted_alignments[j][0] < (i + 1) * line_length]
        space_amounts = [sorted_alignments[index][0] % line_length for index in new_read_indices]
        new_reads = [sorted_reads[index] for index in new_read_indices]
        new_reads_with_spaces = [' ' * space_amounts[j] + new_reads[j] for j in range(len(new_reads))]
        empty_active_read_indices = [index for index in range(len(active_reads)) if active_reads[index] == '']
        for j in range(min(len(new_reads_with_spaces), len(empty_active_read_indices))):
            active_reads[empty_active_read_indices[j]] = new_reads_with_spaces[j]

        if len(new_reads_with_spaces) > len(empty_active_read_indices):
            active_reads += new_reads_with_spaces[len(empty_active_read_indices):]
        printed_reads = ['Read: ' + read[:line_length] for read in active_reads]
        active_reads = [read[line_length:] for read in active_reads]
        while len(active_reads) > 0:
            last_thing = active_reads.pop()
            if last_thing != '':
                active_reads.append(last_thing)
                break
        output_lines = ['Ref:  ' + next_ref] + printed_reads
        output_str += 'Reference index: ' + str(i * line_length) + \
                      '\n' + '\n'.join(output_lines) + '\n\n' + '-' * (line_length + 6) + '\n\n'
    # print output_str
    return output_str


if __name__ == "__main__":
#    folder = 'practice_W_1' #PRACTICE FOLDER

    folder = 'hw1_W_2' #GRADED HOMEWORK
    f_base = '{}_chr_1'.format(folder)

    reads_fn = join(folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)#[:300]
    # This is for speed;
    # If you want to read everything
    # remove the [:300] part of the above line.

    reference_fn = join(folder, 'ref_{}.txt'.format(f_base))
#    reference_fn = join(folder, 'test_ref.txt'.format(f_base)) #TEST REFERENCE FILE
    reference = read_reference(reference_fn)
    # donor_fn = join(folder, 'donor_{}'.format(f_base))
    # donor = read_reference(donor_fn)
    alignments, reads = trivial_algorithm(input_reads, reference)
#    print alignments
#    print reads

    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(folder, 'aligned_{}.txt'.format(f_base))
    # print output_fn

    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
