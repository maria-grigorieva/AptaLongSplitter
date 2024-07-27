import argparse
import os
import uuid
import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import difflib
import numpy as np
from rapidfuzz import fuzz
import re
import logging
from scipy.signal import find_peaks

logger = logging.getLogger()
logging.basicConfig(filename='logfile.log', encoding='utf-8', level=logging.INFO)
formatter = logging.Formatter('%(message)s')
#
# file_handler = logging.FileHandler('logfile.log')
# file_handler.setFormatter(formatter)
# logger.addHandler(file_handler)

parser = argparse.ArgumentParser(
    description='Split long sequencess by barcodes')

# Define the arguments
args_info = [
    ['-i', '--input', str, 'Path to the input fastq file', 'input_data'],
    ['-bar', '--barcodes', str, 'List of barcodes (including complementary)', None]
]

# Add arguments to the parser
for arg_info in args_info:
    parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

# Parse the arguments
args = parser.parse_args()

FILE_PATH = args.input
BARCODES = args.barcodes.split(',')
#
# BARCODE_REGEXP = r'AAGGTTAA.{20,24}CAGCACCT'

def find_fuzzy_substring_matches(rec, barcodes, treshold, result):
    barcode_detected = False
    sequence = str(rec.seq)
    phred_score = rec.letter_annotations["phred_quality"]
    indexes = []
    for id,b in enumerate(barcodes):
        #logger.info(f'Barcode: {b}')
        length = len(b)
        substrings = [sequence[i:i + length] for i in range(len(sequence) - length + 1)]
        arr = [(sequence.find(i), fuzz.ratio(i, b)) for i in substrings if fuzz.ratio(i, b) >= treshold * 100]
        if len(arr) > 0:
            idx, values = zip(*arr)
            peak_indices, _ = find_peaks(values, distance=length)
            peaks = [idx[item] for item in peak_indices.tolist()]
            for p in peaks:
                indexes.append((p, p+len(b)))

        if len(indexes) > 0:
            logger.info(f'Detected {len(indexes)} different barcodes in one line: {indexes}')
            barcode_detected = True
            for i in range(0,len(indexes)):
                start_idx = indexes[i][1]
                end_idx = indexes[i+1][0] if len(indexes) > 1 else len(sequence)+1
                result.append({'sequence': sequence[start_idx:end_idx],
                               'score': phred_score[start_idx:end_idx]})
    if not barcode_detected:
        result.append({'sequence': sequence,
                        'score': phred_score})


def pairwise_similarity_search(rec, barcodes, treshold, result):
    barcode_detected = False
    sequence = str(rec.seq)
    phred_score = rec.letter_annotations["phred_quality"]
    indexes = []
    for id,b in enumerate(barcodes):
        #logger.info(f'Barcode: {b}')
        alignments = pairwise2.align.localms(b, sequence, 5, -4, -4, -.9, one_alignment_only=1)
        logger.info(format_alignment(*alignments[0]))
        if alignments[0].score >= len(b)*5*treshold:
            logger.info(f'Barcode {id} detected')
            indexes.append((alignments[0].start, alignments[0].end))
        if len(indexes) > 0:
            logger.info(f'Detected {len(indexes)} different barcodes in one line: {indexes}')
            barcode_detected = True
            for i in range(0,len(indexes)-1):
                start_idx = indexes[i][1]
                end_idx = indexes[i+1][0] if len(indexes) > 1 else len(sequence)+1
                result.append({'sequence': sequence[start_idx:end_idx],
                               'score': phred_score[start_idx:end_idx]})
    if not barcode_detected:
        result.append({'sequence': sequence,
                        'score': phred_score})


# def barcode_splitter_regexp(rec, barcodes, result):
#     # for b in barcodes:
#     sequence = str(rec.seq)
#     phred_score = rec.letter_annotations["phred_quality"]
#     matches = [(match.start(), match.group()) for match in re.finditer(BARCODE_REGEXP, sequence)]
#     indexes = [(match[0],match[0]+len(match[1])) for match in matches]
#
#     if len(indexes) > 0:
#         new_sequence = sequence[indexes[0][1]:]
#         new_phred_score = phred_score[indexes[0][1]:]
#         result.append({'sequence': new_sequence,
#                        'score': new_phred_score})
#     else:
#         result.append({'sequence': sequence,
#                        'score': phred_score})

def barcode_splitter(rec, barcodes, result):
    barcode_detected = False
    for b in barcodes:
        if not barcode_detected:
            sequence = str(rec.seq)
            phred_score = rec.letter_annotations["phred_quality"]
            substrings = [sequence[i:i+len(b)] for i in range(len(sequence) - len(b) + 1)]
            sub_phred_scores = [phred_score[i:i+len(b)] for i in range(len(phred_score) - len(b) + 1)]
            # averages = [np.median(i) for i in sub_phred_scores]
            # mean_phred = np.mean(averages)
            # remove sequences with low average phred_score value
            # filtered_substrings = []
            # for i in range(0, len(substrings)):
            #     if averages[i] >= mean_phred:
            #         filtered_substrings.append(substrings[i])
            if len(substrings) > 0:

                matches = difflib.get_close_matches(word=b,
                                                    possibilities=substrings,
                                                    n=len(substrings),
                                                    cutoff=0.85)
                if len(matches) > 0:
                    barcode_detected = True
                    # get indexes and sort in an ascending order
                    sorted_matches = sorted([(substring, sequence.find(substring)) for substring in matches], key=lambda x: x[1])

                    sequences_to_remove = [{'sequence': sorted_matches[i][0],
                                           'phred_score': phred_score[sorted_matches[i][1]:sorted_matches[i][1]+len(sorted_matches[i][0])]} for i in range(0,len(sorted_matches))]

                    # initialization
                    s = sorted_matches[0][0]
                    idx = sorted_matches[0][1]
                    ratio = difflib.SequenceMatcher(None, b, s).ratio()

                    X = []

                    # Calculate similarity using SequenceMatcher
                    for item in sorted_matches[1:]:
                        if (item[1] - idx) <= len(b):
                            similarity = difflib.SequenceMatcher(None, b[-8:], item[0][-8:]).ratio()
                            if similarity >= ratio:
                                s = item[0]
                                idx = item[1]
                                ratio = similarity
                        elif (item[1] - idx) > len(b):
                           X.append((str,idx,idx+len(s)))
                           # reinitialization
                           s = item[0]
                           idx = item[1]
                           ratio = difflib.SequenceMatcher(None, b[-8:], s[-8:]).ratio()
                    X.append((s,idx,idx+len(s)))

                    for i in range(0,len(X)):
                        start = X[i][2]
                        end = X[i+1][1] if len(X) > i+1 else len(sequence)
                        result.append({'sequence': sequence[start:end],
                                       'score': phred_score[start:end]})

    if not barcode_detected:
        result.append({'sequence': sequence,
                       'score': phred_score})



def get_all_occurrences(reference, sequence):
    occurrences = fuzzy_substring_search(sequence, reference)
    if len(occurrences)>0:
        return occurrences
    else:
        return -1


def fuzzy_substring_search(s, reference, threshold = 85):
    length = len(reference)
    substrings = [s[i:i + length] for i in range(len(s) - length+1)]
    tmp = [s.find(i) for i in substrings if fuzz.ratio(i, reference) >= threshold]
    return tmp


def reverse_compl_sequence(records):
    new_records = []
    for rec in records:
        occurrences = get_all_occurrences('GGCTTCTGGACTACCTATGC', rec['sequence'])
        if occurrences != -1 and len(occurrences) > 0 and occurrences[0] <= 50:
            new_records.append({'sequence': str(Seq(rec['sequence']).complement())[::-1],
                    'score': rec['score'][::-1]})
        else:
            new_records.append(rec)
    return new_records


def get_sequences():
    """
    Merge sequences from a specified directory
    """
    records = list(SeqIO.parse(FILE_PATH, "fastq"))
    result = []
    counter = 0

    for rec in tqdm.tqdm(records, total=len(list(records)), leave=False, desc=FILE_PATH):
        #barcode_splitter_regexp(rec, BARCODES, result)
        #logger.info('--------------------------------------')
        #logger.info(str(rec.seq))
        find_fuzzy_substring_matches(rec, BARCODES, 0.75, result)

    print(f'Number of detected barcodes: {counter}')

    #result = reverse_compl_sequence(result)

    directory, filename = os.path.split(FILE_PATH)
    filename_base, file_extension = os.path.splitext(filename)
    new_filename = f"{filename_base}_split{file_extension}"

    output_file = os.path.join(directory, new_filename)

    generate_fastq_file(result, output_file)


def generate_fastq_file(data, output_filename):
    records = []
    for i, item in enumerate(data):
        sequence = item['sequence']
        score = item['score']
        record_id = str(uuid.uuid4())
        record = SeqIO.SeqRecord(Seq(sequence), id=record_id, name='', description='')
        record.letter_annotations["phred_quality"] = score
        records.append(record)

    SeqIO.write(records, output_filename, "fastq")
    return output_filename

if __name__ == '__main__':
   get_sequences()

