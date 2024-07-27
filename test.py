import pandas as pd
from scipy.signal import find_peaks
from rapidfuzz import fuzz
import regex
from Bio import SeqIO
from Bio import pairwise2
import time
from Bio.pairwise2 import format_alignment
import numpy as np
from rapidfuzz.process import cdist
import matplotlib
print(matplotlib.rcsetup.all_backends)

matplotlib.use('macosx')
import matplotlib.pyplot as plt

# # Sample DataFrame A
# data = {'index': range(13), 'proportion': [0, 0, 10, 15, 20, 9, 3, 0, 0, 50, 100, 30, 0]}
# df = pd.DataFrame(data)
#
# # Find peaks (local maximums)
# peak_indices, _ = find_peaks(df['proportion'])
#
# # Find valleys (local minimums)
# valley_indices, _ = find_peaks(-df['proportion'])
#
# extremums = [{'index': index, 'proportion': df['proportion'][index]} for index in peak_indices]
# extremums.extend([{'index': index, 'proportion': df['proportion'][index]} for index in valley_indices])
#
# print(extremums)
FILE_PATH = 'input/barcode05_merged.fastq'
s = "CACGCCCGCGGACTATCCCTGTGACAGGAAAAGGTACGGGCGATTTGGCAAACTA"
s_long = "CCGGATTCACTAAGCGGGCGCCGTCCTACGACCCGCGCGCTTTCAGGACCACTCGGGCACGTGGCAGGTCGCTTGCACGCCCGCAGCTGTCCCTGTGACAGGAAAAGGTACGGGCCATTTGGCGAAGCTGGAGCGAGCCTCGGCGGAAAGCTGGGAGGCGCCACCCGGCTTCTGCACAAAGGGCCATCCGGTCGCACAGGGCAGCGGCGCTGCCGGAGGACCAGGGCCGGCGTGCGGCGTCCGGCGAGATCGCGAGCTGCCTCAGGCCGGCGCGCCGCGCTGGGATCCGCCGAGGACCCCGTCGGCGGGAACACCCCGCCAGCTCGCCGAGCTCCGCCGCGCGCCGGCCCCGCCCCGCGCGCTCTAGCTTCGGGTCCTCGGCTCCGCCCCGCTCTGCGACCCCACCACCGCCGTCCCGTGCCCCTGGCCCCGCCCCGCGCCGGATATGCTGGGACGCCCGCGCCTAGAACGCTTTCGTCCCGACGCCCGCAGGTCCTCGCAGTGCGCACCGTTTTGCGACTTGGTGAATTTCTGGTAGCCTCGCTCCCGGAAGAGTGCGAGCTGTCCCTCGGGACGGTGGCGGCCTGGTGGTTCTGCAGGCGCCCTCGCTTCGCCGTCGGTGTGGGCGGCCTGACCCCCACCCATCCGGGCGAGCTTCCGGTGCGCCCAAGTGCCTCCCGGTGTTCCCAGCCTTTCCCGGGCCTGGGGTTGCCTGGACTAGGCTGCGCTGCAGTGACTGTGGACTGGCGTGTGGCGGGGGTGGTGAGGTTAGCCTTACCTCTAGGTGTAGGTGCTGCTTGTCCAGGTTTTGTGTAACCTTTGCTCACAATCATC"
reference = "TGACAGGAAAAGGTACGGGC"
demo_ref = "GGAAAAGGTACGGGC"

def find_fuzzy_substring_matches(s, reference, treshold=75):
    length = len(reference)
    substrings = [s[i:i + length] for i in range(len(s) - length + 1)]
    idx,values = zip(*[(s.find(i),fuzz.ratio(i, reference)) for i in substrings if fuzz.ratio(i, reference) >= treshold])
    peak_indices, _ = find_peaks(values, distance=length)
    return [idx[item] for item in peak_indices.tolist()] if len(peak_indices) > 0 else [-1]

def find_regexp(s, reference, treshold=0.75):
    limit = round(len(reference)*(1-treshold))
    fuzzy_pattern = f'({reference}){{e<={limit}}}'
    match = regex.search(fuzzy_pattern, s, regex.BESTMATCH)
    return [match.span()[0]] if match is not None else [-1]
    #return [s.find(match)] if match is not None else [-1]

    #print(f'{s} and {fuzzy_pattern} = {match}')

def pairwise_similarity_search(s, reference, treshold=0.75):
    alignments = pairwise2.align.localms(reference, s, 5, -4, -4, -.9, one_alignment_only=1)
    if len(alignments) > 0:
        print(format_alignment(*alignments[0]))
        return [alignments[0].start] if alignments[0].score >= len(reference)*5*treshold else [-1]


def get_all_occurrences(reference, all_sequences):
    occurrences = []
    for s in all_sequences:
        current_occurrences = find_fuzzy_substring_matches(s, reference)
        occurrences.extend(current_occurrences)
    unique_values, counts = np.unique(occurrences, return_counts=True)
    data = []
    data_absolute = []
    for value, count in zip(unique_values, counts):
        if value > 0:
            data.append({'index': value, 'proportion': count / len(all_sequences)})
            data_absolute.append({'index': value, 'n_occurrences': count})
    data_df = pd.DataFrame(data)
    print(data_df)
    data_absolute_df = pd.DataFrame(data_absolute)
    print(data_absolute_df)
    #return sorted([i if i>=0 else -1 for i in occurrences])

length = len(reference)
substrings = [s_long[i:i + length] for i in range(len(s_long) - length + 1)]
idx, values = zip(*[(s_long.find(i), fuzz.ratio(i, reference)) for i in substrings])
print(values)

# Plotting the values
plt.plot(values)
plt.xlabel('Position in sequence')
plt.ylabel('Percentage of similarity with the reference')
# plt.title('Plot of the List Values')
plt.show()

# find_regexp(s, reference)

# records = list(SeqIO.parse(FILE_PATH, "fastq"))
# sequences = [str(rec.seq) for rec in records]
#
# start_time = time.time()
# result = cdist([reference], sequences, score_cutoff=0.75)
# print(f"Regexp search for sequences: {(time.time() - start_time)}")
#
# print(result)
#
# # quality test
# LIMITS = [500, 1000, 5000, 10000, 20000, 30000, 40000]
#
# get_all_occurrences(reference, sequences)

# Preformance test
# LIMITS = [500, 1000, 5000, 10000, 20000, 30000, 40000]
#
# for l in LIMITS:
#     records = list(SeqIO.parse(FILE_PATH, "fastq"))[:l]
#     sequences = [str(rec.seq) for rec in records]
#
#     start_time = time.time()
#     for s in sequences:
#         find_regexp(s, reference)
#     print(f"Regexp search for {l} sequences: {(time.time() - start_time)}")
#
#     start_time = time.time()
#     for s in sequences:
#         find_fuzzy_substring_matches(s, reference)
#     print(f"Levenshtein search for {l} sequences: {(time.time() - start_time)}")
#
#     start_time = time.time()
#     for s in sequences:
#         pairwise_similarity_search(s, reference)
#     print(f"BioPython Pairwise2 Similarity search for {l} sequences: {(time.time() - start_time)}")
