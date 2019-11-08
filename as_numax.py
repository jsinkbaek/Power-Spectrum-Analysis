"""
Script used to calculate frequency of maximum oscillations power, by smoothing the power spectrum.
Does not implement the function as_f.smoother, even though it does exactly the same steps.
"""

import matplotlib.pyplot as plt
import as_f
import ps_f
import numpy as np
import math
from pandas import Series

plt.rcParams.update({'font.size': 30})  # Changes font size for all plots

# # Data load from file
# Input data file name
print('Input data to load, exclude extension (.dat)')
# inpt1 = input()  # text input
inpt1 = 'md_sn'
# Call reader function to pull data into list
data = ps_f.dreader(inpt1 + '.dat')

# Preload data lists
freq = []
p = []

# Append data to lists
for row in data:
    freq.append(row[0])
    p.append(row[1])
muHz = 0.000001
freq_muHz = [x / muHz for x in freq]


# # Use whole interval
coord1 = 0
coord2 = len(freq_muHz) - 1
p_select = p[coord1:coord2]

freq_select = freq_muHz[coord1:coord2]
plt.close()

p_array = Series(p_select)

# Set window size in frequency units
wsize_freq = 3 * 25.19  # microHz)

# Find frequency resolution
freq_diff = np.mean(np.diff(freq_select))

# Set window size in index integers
wsize = math.ceil(wsize_freq / freq_diff)

# In case of even window size, set odd window size
if np.fmod(wsize, 2) == 0:
    wsize = wsize - 1

# Smooth by using a hann smoothing window (function call), and convert to list
p_smooth = as_f.hann_smooth(p_array, window_size=wsize).tolist()

# Find maximum values, in order to scale smoothed spectrum for plotting reasons, and to find frequency of max oscill pow
maxp = np.max(p_select)
maxs = np.max(p_smooth)
maxs_indx = np.argmax(p_smooth)

# Calculate scaled spectrum
p_smooth_scaled = [maxp * x / maxs for x in p_smooth]

# Plot smoothed spectrum, normal spectrum, and smoothed+scaled spectrum
plt.figure(1)
plt.plot(freq_select, p_select, 'C0')
plt.plot(freq_select, p_smooth, 'r', linewidth=3)
plt.plot(freq_select, p_smooth_scaled, 'C6--', linewidth=2)
plt.plot(freq_select[maxs_indx], p_smooth[maxs_indx], 'b.', markersize=10)
plt.plot(freq_select[maxs_indx], p_smooth_scaled[maxs_indx], 'b.', markersize=10)
plt.xlabel('Frequency [microHz]')
plt.ylabel('Signal/Noise')
plt.legend(['Signal-to-noise', 'Smoothed', 'Smoothed and scaled', 'Max'])

plt.show(block=True)

print('Max frequency')
print(freq_select[maxs_indx])
