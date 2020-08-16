#!/usr/bin/env python3


import itertools

# this module provides functions for finding the ORI of a DNA sequence

def get_top_repeated_sequences(genome: str, sequence_length: int, number_to_return: int):
    '''
    Runs through all possible permutations of length sequence_length,
    returning the top number_to_return segments from the genome.
    The list returned from the function will contain dictionaries of the following
    format:
    {'sequence': <sequence matched>, 'times_matched': <times matched>}
    '''

    # create dictionary
    sequences_matched = []

    # nucleotide types
    nucleotides = ['A', 'T', 'G', 'C']

    # iterate through all possible combinations of sequence_length
    for sequence in itertools.combinations_with_replacement(nucleotides, sequence_length):
        times_matched = genome.count(''.join(sequence))

        if len(sequences_matched) != 0 and times_matched > sequences_matched[-1]['times_matched']:
            for i in range(len(sequences_matched)-1, -1, -1):
                # iterate backwards through the array, stopping when we find
                # a sequence that's repeated more than the current one
                if sequences_matched[i]['times_matched'] >= times_matched:
                    # splice in the more frequent one, making sure to drop any extras off the end
                    sequences_matched = sequences_matched[:i+1] + \
                                        [{'sequence_matched': sequence, 'times_matched': times_matched}] + \
                                        sequences_matched[i+1:number_to_return-1]
                    break
                elif sequences_matched[i]['times_matched'] < times_matched and i == 0:
                    sequences_matched = [{'sequence_matched': sequence, 'times_matched': times_matched}] + \
                                        sequences_matched[:number_to_return-1]
        elif len(sequences_matched) == 0:
            sequences_matched.append({'sequence_matched': sequence, 'times_matched': times_matched})
    return sequences_matched


if __name__ == '__main__':
    f = open('C:\\Users\\edayley\\Downloads\\Vibrio_cholerae.txt', 'r')
    genome = f.read()
    f.close()
    sequence_length = int(input('sequence length: '))
    for i in get_top_repeated_sequences(genome, sequence_length, 10):
        print(i)