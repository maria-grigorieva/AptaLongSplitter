import pandas as pd
import matplotlib
matplotlib.use("macosx")
import matplotlib.pyplot as plt
from statsmodels.tsa.holtwinters import ExponentialSmoothing
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.signal import find_peaks, peak_prominences, peak_widths, find_peaks_cwt, savgol_filter
import numpy as np
from whittaker_eilers import WhittakerSmoother
from pykalman import KalmanFilter



df = pd.read_csv('result_df.csv')

# def split_list(input_list, chunk_size):
#     to_cut = len(input_list) % chunk_size
#     input_list = input_list[:len(input_list)-to_cut]
#     return [pd.Series(input_list[i:i + chunk_size], index=range(i, i + chunk_size))
#             for i in range(0, len(input_list), chunk_size)]
#
# chunks = split_list(df['proportion'].values, 50)
#
# def get_chunk_extremums(chunks):
#     return [(chunk.idxmax(), chunk.max()) for chunk in chunks]
#
# chunk_max = get_chunk_extremums(chunks)
# print(chunk_max)
# second_values = [t[1] for t in chunk_max]
#
# # Calculate median and maximum
# median_value = np.median(second_values)
# max_value = max(second_values)
#
# print("Median of all second values:", median_value)
# print("Maximum of all second values:", max_value)
#
# pairwise_distances = [abs(chunk_max[i][0] - chunk_max[i+1][0]) for i in range(len(chunk_max) - 1)]
# print(pairwise_distances)
# print(np.mean(pairwise_distances))
# # Create an empty array of zeros with a length of 900
# curve = np.zeros(len(df['proportion'].values))
#
# # Interpolate the given values into the curve
# for tpl in chunk_max:
#     index = tpl[0]
#     value = tpl[1]
#     curve[index] = value
#
# # Perform linear interpolation to fill in the gaps
# indices = np.arange(len(df['proportion'].values))
# curve_interpolated = np.interp(indices, np.where(curve != 0)[0], curve[curve != 0])
#
# # Create pandas series with indexed values of the curve
# series = pd.Series(curve_interpolated, index=indices)
#
# print(series)
# df['smoothed'] = series
# df['smoothed'].plot()
# df['proportion'].plot()
#
#
# peaks, _ = find_peaks(df['smoothed'].values)
# print(peaks)
# plt.scatter(peaks, df.loc[df['index'].isin(peaks), 'proportion'], color='red', zorder=5)
#
# plt.show()
# SMOOTHING
whittaker_smoother = WhittakerSmoother(
    lmbda=len(df), order=2, data_length=len(df)
)

smoothed_data = whittaker_smoother.smooth(df['proportion'])
df['whittaker_smoother'] = smoothed_data
print(df.describe())
df['whittaker_smoother'] = df['whittaker_smoother'].apply(lambda x: 0.0 if x < 0 else x)

df['lowess_smoothed'] = lowess(df['proportion'], range(len(df)), frac=0.06)[:, 1]
# peak_indices, bases = find_peaks(df['proportion'].values, width=1, height=0.0001, distance=20)
# print(peak_indices)
# print(bases)

df['savgol'] = savgol_filter(df['proportion'].values, window_length=40, polyorder=2)

# Define Kalman filter model
kf = KalmanFilter(transition_matrices=[1],
                  observation_matrices=[1],
                  initial_state_mean=df['proportion'].iloc[0],
                  initial_state_covariance=1,
                  observation_covariance=1,
                  transition_covariance=0.01)

# Fit model and make predictions
state_means, _ = kf.filter(df['proportion'])
state_means = state_means.flatten()
df['kalman'] = state_means

peaks, _ = find_peaks(df['whittaker_smoother'].values)
print(peaks)
prominences = peak_prominences(df['whittaker_smoother'].values, peaks)[0]
results_full = peak_widths(df['whittaker_smoother'].values, peaks, rel_height=1)[0]
avg_prominence = np.median(prominences)
print('PROMINENCES')
print(prominences)
print(avg_prominence)
print('WIDTHS')
print(results_full)
avg_width = np.percentile(results_full, 15)
print(avg_width)
peak_indices_smoothed, bases_smoothed = find_peaks(df['whittaker_smoother'].values,
                                                   #wlen=20,
                                                   width=avg_width,
                                                   # height=0.0001,
                                                   distance=20,
                                                   prominence=avg_prominence)
                                                   #prominence=0.0001)
print(peak_indices_smoothed)
print(bases_smoothed)

# df['savgol'] = savgol_filter(df['proportion'], window_length=10, polyorder=1, mode='nearest')
# peak_indices_savgol, bases_savgol = find_peaks(df['savgol'].values, width=1, height=0.0001,distance=20)
# print(peak_indices_savgol)
# print(bases_savgol)
#
# peak_indices_ = find_peaks_cwt(df['smoothed'].values, np.arange(1,200))
# print(peak_indices_)
df['proportion'].plot()
# df['smoothed'].plot()
df['whittaker_smoother'].plot()
plt.scatter(peak_indices_smoothed, df.loc[df['index'].isin(peak_indices_smoothed), 'proportion'], color='red', zorder=5)

# df['smoothed'].plot()
# df['savgol'].plot()

# df['proportion'].ewm(span=50).mean().plot()

plt.show()
