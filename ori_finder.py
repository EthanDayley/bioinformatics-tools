#!/usr/bin/env python3


from typing import Dict, List

# this module provides functions for finding the ORI of a DNA sequence

def get_top_repeated_kmers(genome: str, kmer_length: int, max_num_matches: int):
    '''
    Finds the top repeated k-mers in the genome of length kmer_length,
    returning a list of length <= max_num_matches
    '''

    # initialize dict for storing k-mers
    kmers: Dict[str,int] = {}

    # for efficiency, variable stores maximum index to loop to
    max_box_start = len(genome) - kmer_length

    # loop through genome, sliding a "box" of length kmer_length along it
    for box_start in range(max_box_start):
        curr_kmer = genome[box_start:box_start+kmer_length]
        if curr_kmer in kmers:
            kmers[curr_kmer] += 1
        else:
            kmers[curr_kmer] = 1

    # initialize list of top matches
    top_kmers: List[Dict] = []

    # loop through kmers, creating sorted list of top matches
    for kmer in kmers:
        if len(top_kmers) == 0 or \
            top_kmers[-1]['num_occurrences'] >= kmers[kmer] and len(top_kmers) < max_num_matches:
            # list empty OR list below target length and num_ocurrences < last element
            # NOTE: case doesn't capture equal num_ocurrences if length >= max_num_matches
            top_kmers.append({'kmer': kmer, 'num_occurrences': kmers[kmer]})
        elif top_kmers[-1]['num_occurrences'] <=  kmers[kmer]:
            # loop backwards through top_kmers
            for i in range(len(top_kmers)-1, -1, -1):
                if top_kmers[i]['num_occurrences'] >= kmers[kmer]:
                    # define new index (for readability)
                    new_index = i + 1

                    # break if smaller than last element
                    if new_index >= max_num_matches:
                        break

                    # insert new value into list
                    top_kmers = top_kmers[:new_index] \
                        + [{'kmer': kmer, 'num_occurrences': kmers[kmer]}] \
                        + top_kmers[new_index:max_num_matches]
                    break

                elif top_kmers[i]['num_occurrences'] < kmers[kmer] and i == 0:
                    # insert new value into list
                    top_kmers = [{'kmer': kmer, 'num_occurrences': kmers[kmer]}] \
                                + top_kmers[:max_num_matches]
                    break

    return top_kmers



if __name__ == '__main__':
    f = open('C:\\Users\\edayley\\Downloads\\Vibrio_cholerae.txt', 'r')
    genome = f.read()
    f.close()
    sequence_length = int(input('sequence length: '))
    for i in get_top_repeated_kmers(genome, sequence_length, 50):
        print(i)