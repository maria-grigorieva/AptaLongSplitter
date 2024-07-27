import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from rapidfuzz import fuzz
import numpy as np
import plotly.graph_objects as go
import re
from tqdm import tqdm
from scipy.signal import find_peaks, peak_prominences, peak_widths
import json
import logging
import configparser
import os
import itertools
import regex
from statsmodels.nonparametric.smoothers_lowess import lowess
from whittaker_eilers import WhittakerSmoother

# Read parameters from the config file
config = configparser.ConfigParser()
config.read('config.ini')

FILE_PATH = config['Parameters']['INPUT_FILE']
filename = os.path.splitext(os.path.basename(FILE_PATH))[0]
filetype = os.path.splitext(os.path.basename(FILE_PATH))[1]
LOG_PATH = os.path.join(os.path.dirname(FILE_PATH), filename + ".log")
JSON_PATH = os.path.join(os.path.dirname(FILE_PATH), filename + ".json")

FUZZY = bool(config['Parameters']['FUZZY'])
TRESHOLD = float(config['Parameters']['TRESHOLD'])
LIMIT = int(config['Parameters']['LIMIT'])
SMOOTHING = config['Parameters']['SMOOTHING']

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)  # Set the logging level

formatter = logging.Formatter('%(message)s')  # Define the formatter

file_handler = logging.FileHandler(LOG_PATH, encoding='utf-8')  # Create a FileHandler
file_handler.setFormatter(formatter)  # Set the formatter for the FileHandler
logger.addHandler(file_handler)  # Add the FileHandler to the logger


logger.info(FILE_PATH)
logger.info(LOG_PATH)
logger.info(JSON_PATH)
logger.info(FUZZY)
logger.info(TRESHOLD)


SEQUENCES = []
if filetype == '.fastq':
    for key, value in config['Sequences'].items():
        param = {'type': key, 'sequence': value, 'occurences': []}
        SEQUENCES.append(param)

        if config['Parameters']['COMPLEMENTARY'] == 'True':
            seq_obj = Seq(value)
            reverse_complement = seq_obj.reverse_complement()
            reverse_complement_str = str(reverse_complement)
            param_revcompl = {'type': f'{key}_REVCOMPL', 'sequence': reverse_complement_str, 'occurences': []}
            SEQUENCES.append(param_revcompl)

elif filetype == '.json':
    f = open(FILE_PATH)
    PRIMERS = json.load(f)

def pairwise_similarity_search(s, reference, treshold):
    alignments = pairwise2.align.localms(reference, s, 5, -4, -4, -.9, one_alignment_only=1)
    if len(alignments) > 0:
        logger.info(format_alignment(*alignments[0]))
        return [alignments[0].start] if alignments[0].score >= len(reference)*5*treshold else []
    else:
        return [-1]


def find_fuzzy_substring_matches(s, reference, treshold):
    length = len(reference)
    substrings = [s[i:i + length] for i in range(len(s) - length + 1)]
    arr = [(s.find(i),fuzz.ratio(i, reference)) for i in substrings if fuzz.ratio(i, reference) >= treshold*100]
    if len(arr) > 0:
        idx, values = zip(*arr)
        peak_indices, _ = find_peaks(values, distance=length)
        return [idx[item] for item in peak_indices.tolist()]
    else:
        return [-1]

def find_fuzzy_regex(s, reference, treshold):
    limit = round(len(reference) * (1 - treshold))
    fuzzy_pattern = f'({reference}){{e<={limit}}}'
    match = regex.search(fuzzy_pattern, s, regex.BESTMATCH)
    return [match.span()[0]] if match is not None else [-1]

def get_all_occurrences(reference, type, all_sequences):
    occurrences = []
    logger.info(f'Searching for {reference} in all sequences...')
    for s in tqdm(all_sequences, desc=f"Searching for {type}: {reference}", unit="sequence"):
        # current_occurrences = find_fuzzy_regex(s, reference, TRESHOLD)
        #current_occurrences = pairwise_similarity_search(s, reference, TRESHOLD)
        current_occurrences = find_fuzzy_substring_matches(s, reference, TRESHOLD) if FUZZY \
                            else [m.start() for m in re.finditer(reference, s)]
        occurrences.extend(current_occurrences if len(current_occurrences) > 0 else [-1])
    return sorted(occurrences)

def moving_average(data, window_size):
    return data.rolling(window=window_size).mean()

def smooth_data(df, smooth_type='whittaker'):
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
    return target


def get_peak_occurrences(x):
    unique_values, counts = np.unique(x['occurences'], return_counts=True)
    data, data_absolute = [], []

    for value, count in zip(unique_values, counts):
        if value > 0:
            data.append({'index': value, 'proportion': count / len(sequences)})
            data_absolute.append({'index': value, 'n_occurrences': count})
    if len(data) > 0:
        df = pd.DataFrame(data)
        df_abs = pd.DataFrame(data_absolute)
        all_indexes = pd.Series(range(0, avg_length))
        result_df = all_indexes.to_frame('index').merge(df, on='index', how='left').fillna(0)
        result_df.to_csv('result_df.csv')
        result_df_abs = all_indexes.to_frame('index').merge(df_abs, on='index', how='left').fillna(0)

        target = smooth_data(result_df, SMOOTHING)

        peaks, initial_bases = find_peaks(result_df[target].values)
        prominences = peak_prominences(result_df[target].values, peaks)[0]
        avg_prominence = np.mean(prominences)
        widths = peak_widths(result_df[target].values, peaks, rel_height=1)[0]
        peak_indices, bases = find_peaks(result_df[target].values,
                                         width=np.percentile(widths,15),
                                         # distance=round(len(x['sequence'])),
                                         # height=0.0001,
                                         prominence=avg_prominence)
        if len(peak_indices) == 0:
            peak_indices = peaks
            bases = initial_bases
        extremums = []
        for i in range(0, len(peak_indices)):
            peak_index = peak_indices[i]
            extremums.append(aggregate_peak_values(i, result_df, result_df_abs, peak_index, bases))

        x['peaks'] = extremums if len(extremums) > 0 else []
        x['value_counts'] = result_df.to_dict('records') if len(data) > 0 else []
        x['value_counts_abs'] = result_df_abs.to_dict('records') if len(data) > 0 else []


def aggregate_peak_values(step, result_df, result_df_abs, peak_index, bases):
    peak_proportion = result_df['proportion'][peak_index]
    peak_occurrences = result_df_abs['n_occurrences'].iloc[peak_index]
    left_bases = bases['left_bases'][step]
    right_bases = bases['right_bases'][step]
    total_proportion = np.round(np.sum(result_df.iloc[left_bases:right_bases]['proportion'].values), 4)
    total_occurrences = np.round(np.sum(result_df_abs.iloc[left_bases:right_bases]['n_occurrences'].values), 4)

    return {'peak_index': peak_index,
          'peak_proportion': peak_proportion,
          'peak_occurrences': peak_occurrences,
          'left_bases': left_bases,
          'right_bases': right_bases,
          'total_proportion': total_proportion,
          'total_occurrences': total_occurrences}

def plot_distribution_proportions(limit):
    fig = go.Figure()
    for p in SEQUENCES:
        title = p['type']
        if 'value_counts' in p:
            df = pd.DataFrame(p['value_counts']).head(limit)
            total_proportion = np.round(np.sum([item['total_proportion'] for item in p['peaks']]),2) if len(p['peaks']) > 0 else 0

            peaks_info = ""
            for i,item in enumerate(p['peaks']):
                idx = item['peak_index']
                proportion = item['total_proportion']
                peaks_info += f'{idx}({proportion})  '
            fig.add_trace(go.Scatter(
                x=df['index'],
                y=df['proportion'],
                mode='lines',
                name=f'{title}: Total Proportion = {total_proportion}, PEAKS: {peaks_info}',
                )
            )
            if len(SEQUENCES) == 1:
                # Add dots at specific x positions
                peak_positions = [i['peak_index'] for i in p['peaks']]
                fig.add_trace(go.Scatter(
                    x=peak_positions,
                    y=[df.loc[df['index'] == x, 'proportion'].values[0] for x in peak_positions],
                    text=[i['total_proportion'] for i in p['peaks']],
                    textposition='top center',
                    textfont=dict(color='red'),
                    mode='markers+text',
                    marker=dict(size=8, color='red'),
                    showlegend=False
                ))
                # Smoothing line
                x_smooth = df['index']
                y_smooth = df['smoothed']
                fig.add_trace(go.Scatter(
                    x=x_smooth,
                    y=y_smooth,
                    mode='lines',
                    showlegend=False
                ))

    fig.update_layout(
        width=1000,
        height=600,
        legend=dict(
            x=0,
            y=1.0,
            xanchor='left',
            yanchor='bottom'
        ),
        barmode='overlay',
        title={
            'text': "Positional distribution of occurrences proportional to the total number of sequences",
            'y': 0.05,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'bottom'}
    )
    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)

    return fig

def plot_distribution_absolute(limit):
    fig = go.Figure()
    for p in SEQUENCES:
        title = p['type']
        if 'value_counts_abs' in p:
            df = pd.DataFrame(p['value_counts_abs']).head(limit)
            n_occurrences = np.round(np.sum([item['total_occurrences'] for item in p['peaks']]),2) if len(p['peaks']) > 0 else 0
            peaks_info = ""
            for item in p['peaks']:
                idx = item['peak_index']
                occurrences = item['total_occurrences']
                peaks_info += f'{idx}({occurrences})  '
            fig.add_trace(go.Scatter(
                x=df['index'],
                y=df['n_occurrences'],
                mode='lines',
                name=f'{title}: Total Amount = {n_occurrences}, PEAKS: {peaks_info}',
            )
            )
    fig.update_layout(
        width=1000,
        height=600,
        legend=dict(
            x=0,
            y=1.0,
            xanchor='left',
            yanchor='bottom'
        ),
        barmode='overlay',
        title={
            'text': "Positional distribution of occurrences in absolute numbers",
            'y': 0.05,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'bottom'}
    )
    fig.update_layout(barmode='overlay')
    if len(fig.data) == 1:
        fig.update_traces(showlegend=True, name=f'{title}: Total Amount = {n_occurrences}, PEAKS: {peaks_info}',)
    fig.update_traces(opacity=0.75)
    return fig


if filetype == '.fastq':

    records = list(SeqIO.parse(FILE_PATH, "fastq"))
    print(len(records))

    sequences = [str(rec.seq) for rec in records][:LIMIT if LIMIT>0 else len(records)]
    avg_length = int(np.mean([len(s) for s in sequences]))
    logger.info(f'Number of sequences {len(sequences)}')

    # BARCODE['occurences'] = get_all_occurrences(BARCODE['sequence'], BARCODE['type'], sequences)
    # get_peak_occurrences(BARCODE)

    for p in SEQUENCES:
        p['occurences'] = get_all_occurrences(p['sequence'], p['type'], sequences)
        get_peak_occurrences(p)

    max_right_bases = max(
        int(p['right_bases']) for p in itertools.chain.from_iterable(d['peaks'] for d in SEQUENCES)) + 50

    json.dump(SEQUENCES, open(JSON_PATH, 'w'), default=str)

    # PRIMERS.append(BARCODE)
    fig1 = plot_distribution_proportions(max_right_bases)
    fig2 = plot_distribution_absolute(max_right_bases)
    fig1.show()
    fig2.show()

elif filetype == '.json':
    max_right_bases = max(
        int(p['right_bases']) for p in itertools.chain.from_iterable(d['peaks'] for d in SEQUENCES))

    fig1 = plot_distribution_proportions(max_right_bases)
    fig2 = plot_distribution_absolute(max_right_bases)
    fig1.show()
    fig2.show()
