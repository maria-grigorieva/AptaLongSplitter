import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import Align
from Bio.motifs import Motif

from rapidfuzz import fuzz
import numpy as np
from tqdm import tqdm
from scipy.signal import find_peaks, peak_prominences, peak_widths, savgol_filter
import configparser
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
from whittaker_eilers import WhittakerSmoother
from statsmodels.nonparametric.kernel_regression import KernelReg
from confsmooth import confsmooth
from scipy.stats import variation
from Bio.Align import AlignInfo


def find_fuzzy_substring_matches(s, reference, treshold):
    length = len(reference)
    substrings = [s[i:i + length] for i in range(len(s) - length + 1)]
    arr = [(s.find(i),fuzz.ratio(i, reference)) for i in substrings if fuzz.ratio(i, reference) >= treshold*100]
    if len(arr) > 0:
        idx, values = zip(*arr)
        peak_indices, _ = find_peaks(values, distance=length)
        if len(peak_indices) > 0:
            # NEW!!!
            return [(idx[item], s[idx[item]:idx[item]+length]) for item in peak_indices.tolist()]
        else:
            return None
    else:
        return None

def get_all_occurrences(reference, type, all_sequences, n_records, avg_length, threshold=0.75):
    occurrences = []
    for s in tqdm(all_sequences, desc=f"Searching for {type}: {reference}", unit="sequence"):
        try:
            #current_occurrences = find_fuzzy_substring_matches(s, reference, threshold)
            occurrences.extend(find_fuzzy_substring_matches(s, reference, threshold))
        except Exception as e:
            pass
        # occurrences.extend(current_occurrences if len(current_occurrences) > 0 else [-1])
    if len(occurrences) > 0:
        # NEW !!!!
        arr = np.array(sorted(occurrences), dtype=[('index', int), ('reference', 'O')])
        arr_df = pd.DataFrame(arr)
        unique_values, counts = np.unique(arr_df['index'], return_counts=True)
        # unique_counts = arr_df['index'].value_counts()
        consensus_values = []
        for unique_occurrence in unique_values:
            subset = arr_df[arr_df['index'] == unique_occurrence]['reference'].values
            records = [SeqRecord(Seq(seq), id=f"seq{i + 1}") for i, seq in enumerate(subset)]
            msa = MultipleSeqAlignment(records)
            alignment = msa.alignment
            motif = Motif("ACGT", alignment)
            consensus_values.append(str(motif.consensus))
            # summary_align = AlignInfo.SummaryInfo(msa)
            # consensus = summary_align.dumb_consensus()
            # consensus_values.append(str(consensus))

        data = [{'index': value, 'reads': count, 'proportion': round(count / n_records, 4), 'consensus': consensus} for value, count, consensus in
                zip(unique_values, counts, consensus_values)]
        # OLD!!!!
        #unique_values, counts = np.unique(sorted(occurrences), return_counts=True)
        #data = [{'index': value, 'reads': count, 'proportion': round(count / n_records,4)} for value, count in zip(unique_values, counts)]
        df = pd.DataFrame(data)
        all_indexes = pd.Series(range(0, avg_length))
        result = all_indexes.to_frame('index').merge(df, on='index', how='left').fillna(0)
        return result
    else:
        return []


def moving_average(arr, window_size):
    return np.convolve(arr, np.ones(window_size) / window_size, mode='valid')

def calculate_noise(arr, window_size=3):
    smoothed = moving_average(arr, window_size)
    extended_smoothed = np.concatenate((arr[:window_size - 1], smoothed))
    noise = arr - extended_smoothed
    return noise

def calculate_snr(arr, window_size=3):
    noise = calculate_noise(arr, window_size)
    noise_power = np.var(noise)
    signal_power = np.var(arr)

    snr = 10 * np.log10(signal_power / noise_power)
    return snr

def signaltonoise(a, axis=0, ddof=0):
    a = np.asanyarray(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    return np.where(sd == 0, 0, m/sd)

def noise_level(data, reference):
    return variation(data['reads'], axis=0)
    # savgol = savgol_filter(data['reads'], len(reference), 1)
    # return round(np.std(data['reads'] - savgol),4)

def smooth_data(df, reference, smooth_type='whittaker'):
    target = 'proportion'
    if smooth_type.lower() == 'none':
        pass
    else:
        target = 'smoothed'
        if smooth_type.lower() == 'whittaker':
            whittaker_smoother = WhittakerSmoother(
                lmbda=len(df), order=2, data_length=len(df)
            )
            df[target] = whittaker_smoother.smooth(df['proportion'])
        elif smooth_type.lower() == 'lowess':
            df[target] = lowess(df['proportion'], range(len(df)), frac=0.1)[:, 1]
        elif smooth_type.lower() == 'savgol':
            df[target] = savgol_filter(df['proportion'].values, window_length=len(reference), polyorder=2)
        elif smooth_type.lower() == 'confsmooth':
            print('Estimated noise level =', noise_level(df['proportion'], reference))
            df[target] = confsmooth(df['proportion'], noise_level, confidence=0.9, deg=2)
    return target

def get_peak_occurrences(x, smoothing):
    target = smooth_data(x['occurences'], x['sequence'], smoothing)
    peaks, initial_bases = find_peaks(x['occurences'][target].values)
    prominences = peak_prominences(x['occurences'][target].values, peaks)[0]
    avg_prominence = np.mean(prominences)
    widths = peak_widths(x['occurences'][target].values, peaks, rel_height=1)[0]
    peak_indices, bases = find_peaks(x['occurences'][target].values,
                                     #width=np.percentile(widths,15),
                                     # distance=round(len(x['sequence'])),
                                     # height=0.0001,
                                     prominence=avg_prominence)
    if len(peak_indices) == 0:
        peak_indices = peaks
        bases = initial_bases

    if bases:
        extremums = []
        for i in range(0, len(peak_indices)):
            peak_index = peak_indices[i]
            extremums.append(aggregate_peak_values(i,
                                                   x['occurences'],
                                                   peak_index,
                                                   bases))
        # Calculate peak distances
        for i in range(1, len(extremums)):
            extremums[i]['peak_dist'] = extremums[i]['peak_index'] - extremums[i - 1]['peak_index']

        x['peaks'] = extremums if len(extremums) > 0 else []
        x['average_peaks_distance'] = calculate_average_peaks_distance(x['peaks'])
    x['value_counts'] = x['occurences'].to_dict('records') if x['occurences'].shape[0] > 0 else []
    # x['value_counts_abs'] = x['occurences'].to_dict('records') if x['occurences'].shape[0] > 0 else []

def aggregate_peak_values(step, df, peak_index, bases):
    left_bases = bases['left_bases'][step] if 'left_bases' in bases else 0
    right_bases = bases['right_bases'][step] if 'right_bases' in bases else 0
    total_proportion = np.round(np.sum(df.iloc[left_bases:right_bases]['proportion'].values), 4)
    total_occurrences = np.round(np.sum(df.iloc[left_bases:right_bases]['reads'].values), 4)
    consensus = df.loc[int(peak_index),'consensus']
    return {'peak_index': int(peak_index),
          'left_bases': int(left_bases),
          'right_bases': int(right_bases),
          'total_proportion': float(total_proportion),
          'total_reads': int(total_occurrences),
          'motif_consensus': str(consensus)}

def calculate_average_peaks_distance(peaks):
    indexes = [p['peak_index'] for p in peaks]
    if len(indexes) > 2:
        # Initialize a list to hold the distances
        distances = []

        # Calculate distances between consecutive elements
        for i in range(len(indexes) - 1):
            distance = indexes[i + 1] - indexes[i]
            distances.append(distance)

        # Calculate the average distance
        average_distance = sum(distances) / len(distances)
        return average_distance
    else:
        return 0


def main(SESSION):
    # Read parameters from the config file
    config = configparser.ConfigParser()
    config.read(os.path.join(SESSION, 'config.ini'))
    FILE_PATH = config['Parameters']['input_file']
    filename = os.path.splitext(os.path.basename(FILE_PATH))[0]
    filetype = os.path.splitext(os.path.basename(FILE_PATH))[1]
    THRESHOLD = float(config['Parameters']['threshold'])
    LIMIT = int(config['Parameters']['limit'])
    SMOOTHING = config['Parameters']['smoothing']
    # JSON_PATH = os.path.join(FILE_PATH, filename+'.json')
    result_data = {}

    SEQUENCES = []
    for key, value in config['Sequences'].items():
        param = {'type': key, 'sequence': value, 'occurences': []}
        SEQUENCES.append(param)

    records = list(SeqIO.parse(FILE_PATH, "fastq"))
    n_records = len(records)

    sequences = [str(rec.seq) for rec in records][:LIMIT if LIMIT > 0 else n_records]
    avg_length = int(np.mean([len(s) for s in sequences]))

    for p in SEQUENCES:
        p['occurences'] = get_all_occurrences(p['sequence'], p['type'], sequences, n_records, avg_length, THRESHOLD)
        if len(p['occurences']) > 0:
            #p['noise_level'] = calculate_snr(p['occurences']['reads'], round(len(p['sequence'])/2))
            p['noise_level'] = float(np.round(signaltonoise(p['occurences']['reads']),4))
            # p['noise_level'] = noise_level(p['occurences'], p['sequence'])
            p['total_reads'] = int(np.sum(p['occurences']['reads']))
            p['total_proportion'] = float(np.round(np.sum(p['occurences']['proportion']),4))
            get_peak_occurrences(p, smoothing=SMOOTHING)

    result_data['sequences'] = SEQUENCES
    result_data['parameters'] = {'n_records': n_records,
                                 'smoothing': SMOOTHING,
                                 'limit': LIMIT,
                                 'threshold': THRESHOLD,
                                 'file_path': FILE_PATH,
                                 'avg_noise_level': np.round(np.mean([item['noise_level'] for item in SEQUENCES if 'noise_level' in item]),4)
                                 }

    return result_data

