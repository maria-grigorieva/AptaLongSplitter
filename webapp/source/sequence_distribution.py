import pandas as pd
from Bio import SeqIO
from rapidfuzz import fuzz
import numpy as np
from tqdm import tqdm
from scipy.signal import find_peaks, peak_prominences, peak_widths, savgol_filter
import configparser
import os
from statsmodels.nonparametric.smoothers_lowess import lowess
from whittaker_eilers import WhittakerSmoother

def find_fuzzy_substring_matches(s, reference, treshold):
    length = len(reference)
    substrings = [s[i:i + length] for i in range(len(s) - length + 1)]
    arr = [(s.find(i),fuzz.ratio(i, reference)) for i in substrings if fuzz.ratio(i, reference) >= treshold*100]
    if len(arr) > 0:
        idx, values = zip(*arr)
        peak_indices, _ = find_peaks(values, distance=length)
        if len(peak_indices) > 0:
            return [[idx[item] for item in peak_indices.tolist()][0]]
        else:
            return [-1]
    else:
        return [-1]

def get_all_occurrences(reference, type, all_sequences, threshold=0.75):
    occurrences = []
    for s in tqdm(all_sequences, desc=f"Searching for {type}: {reference}", unit="sequence"):
        current_occurrences = find_fuzzy_substring_matches(s, reference, threshold)
        occurrences.extend(current_occurrences if len(current_occurrences) > 0 else [-1])
    return sorted(occurrences)

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
    return target

def get_peak_occurrences(x, sequences, avg_length, smoothing):
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
        # result_df.to_csv('result_df.csv')
        result_df_abs = all_indexes.to_frame('index').merge(df_abs, on='index', how='left').fillna(0)

        target = smooth_data(result_df, x['sequence'], smoothing)

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
            extremums.append(aggregate_peak_values(i,
                                                   result_df,
                                                   result_df_abs,
                                                   peak_index,
                                                   bases))
        # Calculate peak distances
        for i in range(1, len(extremums)):
            extremums[i]['peak_dist'] = extremums[i]['peak_index'] - extremums[i - 1]['peak_index']

        x['peaks'] = extremums if len(extremums) > 0 else []
        x['average_peaks_distance'] = calculate_average_peaks_distance(x['peaks'])
        # x['peaks'][x['type']] = extremums if len(extremums) > 0 else []
        x['value_counts'] = result_df.to_dict('records') if len(data) > 0 else []
        x['value_counts_abs'] = result_df_abs.to_dict('records') if len(data) > 0 else []

def aggregate_peak_values(step, result_df, result_df_abs, peak_index, bases):
    # peak_proportion = result_df['proportion'][peak_index]
    # peak_occurrences = result_df_abs['n_occurrences'].iloc[peak_index]
    left_bases = bases['left_bases'][step]
    right_bases = bases['right_bases'][step]
    total_proportion = np.round(np.sum(result_df.iloc[left_bases:right_bases]['proportion'].values), 4)
    total_occurrences = np.round(np.sum(result_df_abs.iloc[left_bases:right_bases]['n_occurrences'].values), 4)
    return {'peak_index': peak_index,
          'left_bases': left_bases,
          'right_bases': right_bases,
          'total_proportion': total_proportion,
          'total_occurrences': total_occurrences}

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

    SEQUENCES = []
    for key, value in config['Sequences'].items():
        param = {'type': key, 'sequence': value, 'occurences': []}
        SEQUENCES.append(param)

    records = list(SeqIO.parse(FILE_PATH, "fastq"))
    n_records = len(records)
    sequences = [str(rec.seq) for rec in records][:LIMIT if LIMIT > 0 else n_records]
    avg_length = int(np.mean([len(s) for s in sequences]))

    for p in SEQUENCES:
        p['occurences'] = get_all_occurrences(p['sequence'], p['type'], sequences, THRESHOLD)
        get_peak_occurrences(p, sequences=sequences, avg_length=avg_length, smoothing=SMOOTHING)

    # max_right_bases = max(
    #     int(p['right_bases']) for p in itertools.chain.from_iterable(d['peaks'] for d in SEQUENCES)) + 50

    return SEQUENCES
    #json.dump(SEQUENCES, open(JSON_PATH, 'w'), default=str)

    # fig1 = plot_distribution_proportions(SEQUENCES, max_right_bases, SMOOTHING)
    # graphJSON1 = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)
    #
    # peaks_table = make_peaks_subplots(SEQUENCES)
    # peaksJSON = json.dumps(peaks_table, cls=plotly.utils.PlotlyJSONEncoder)
    # # fig2 = plot_distribution_absolute(SEQUENCES, max_right_bases)
    # # graphJSON2 = json.dumps(fig2, cls=plotly.utils.PlotlyJSONEncoder)
    # # detected_peaks = [p['peaks'] for p in SEQUENCES]
    # # # Select specific fields for each dictionary in the list
    # # selected_peaks = []
    # # for peak_list in detected_peaks:
    # #     selected_peak_list = []
    # #     for peak_dict in peak_list:
    # #         selected_peak_dict = {key: peak_dict[key] for key in
    # #                               ['peak_index', 'left_bases', 'right_bases', 'total_proportion',
    # #                                'total_occurrences']}
    # #         selected_peak_list.append(selected_peak_dict)
    # #     selected_peaks.append(selected_peak_list)
    #
    # return graphJSON1, peaksJSON
