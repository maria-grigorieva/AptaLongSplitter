import argparse
import os
import uuid
import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
import difflib

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

def barcode_splitter(rec, barcodes, result):
    barcode_detected = False
    for b in barcodes:
        if not barcode_detected:
            sequence = str(rec.seq)
            phred_score = rec.letter_annotations["phred_quality"]
            substrings = [sequence[i:i+len(b)] for i in range(len(sequence) - len(b) + 1)]
            matches = difflib.get_close_matches(word=b,
                                                possibilities=substrings,
                                                n=len(substrings),
                                                cutoff=0.85)
            if len(matches) > 0:
                barcode_detected = True
                # get indexes and sort in an ascending order
                sorted_matches = sorted([(substring, sequence.find(substring)) for substring in matches], key=lambda x: x[1])

                # initialization
                s = sorted_matches[0][0]
                idx = sorted_matches[0][1]
                ratio = difflib.SequenceMatcher(None, b[-8:], s[-8:]).ratio()

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


def get_sequences():
    """
    Merge sequences from a specified directory
    """
    records = list(SeqIO.parse(FILE_PATH, "fastq"))
    result = []

    for rec in tqdm.tqdm(records, total=len(list(records)), leave=False, desc=FILE_PATH):
        barcode_splitter(rec, BARCODES, result)

    directory, filename = os.path.split(FILE_PATH)
    filename_base, file_extension = os.path.splitext(filename)
    new_filename = f"{filename_base}_splitted{file_extension}"

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

