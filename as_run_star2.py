"""
Script to run through all the power-spectrum analysis for TIC 71134596. Similar to the one for nu Indi, but slightly
modified.
-   Calculates median noise level for 3 intervals selected by graphical input.
-   Does peak-finding for the 3 intervals by using the function detect_peaks (https://github.com/demotu/BMC).
-   Finds peaks above power-noise limit of 6, by using as_f.good_peak_finder.
-   Forces manual assessment of each peak, in order to determine which peak is a primary peak for an oscillation, and
    which is simply a very tall secondary peak (sinc-like behaviour). Result is saved to a .dat-file, and can be
    reloaded instead of repeating the assessment. This overwrites the previous calculations if done.
-   Calculates autocorrelation using only data from peaks in power-spectrum,
    in order to (using graphical input selection) determine large frequency separation.
-   Plots echelle heatmap diagrams multiple times using only data from peaks in power-spectrum, in order to determine
    small frequency separation.
-   Plots echelle heatmap diagrams using all peaks that have or have not gone through manual assessment, in order to
    compare.
"""

import math
import ps_f
import matplotlib.pyplot as plt
import numpy as np
from detect_peaks import detect_peaks
import as_f

plt.rcParams.update({'font.size': 25})  # Changes font size for all plots

# # Data load from file
# Input data file name
print('Input data to load, exclude extension (.dat)')
#inpt1 = input()
inpt1 = 'star2_ps'
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


# # Find data medians (noise level) for all intervals
number_of_intervals = 3  # Sets the amount of medians necessary for the data (intervals)
[median_noise, intv, p_intv, freq_intv] = as_f.median_noise(freq, p, number_of_intervals=number_of_intervals)

# # Pretty plot
pn_limit = 6

plt.subplot(2, 2, 1)
plt.plot(freq_muHz, p, 'C0')
colors = ['C1', 'C2', 'C3']
for i in range(0, len(intv)):
    freq_i = [x / muHz for x in freq_intv[i]]
    plt.plot(freq_i, [median_noise[i]] * len(freq_intv[i]), colors[i])
plt.yscale('log', basey=10)
plt.xscale('log', basex=10)
plt.ylabel('Power Spectrum')

plt.subplot(2, 2, 2)
plt.plot(freq_muHz, p, 'C0')
for i in range(0, len(intv)):
    freq_i = [x / muHz for x in freq_intv[i]]
    plt.plot(freq_i, [median_noise[i]] * len(freq_intv[i]), colors[i])

plt.subplot(2, 2, 3)
p_o_n = as_f.func_over_noise(median_noise, p_intv)
p_o_n_freq = []
for i in range(0, number_of_intervals):
    freq_i = [x / muHz for x in freq_intv[i]]
    p_o_n_freq.extend(freq_i)
plt.plot(p_o_n_freq, p_o_n, 'C0')
plt.plot(p_o_n_freq, [pn_limit] * len(p_o_n_freq), 'k')
plt.yscale('log', basey=10)
plt.xscale('log', basex=10)
plt.ylabel('Power/Noise')
plt.xlabel('Frequency [microHz]')

plt.subplot(2, 2, 4)
plt.plot(p_o_n_freq, p_o_n, 'C0')
plt.plot(p_o_n_freq, [pn_limit] * len(p_o_n_freq), 'k')
plt.xlabel('Frequency [microHz]')

plt.show()

# # Peak finding (peak larger than median)

# Initialize lists
p_peaks = [[] for i in range(0, len(p_intv))]
freq_peaks = [[] for i in range(0, len(freq_intv))]
freq_muHz_peaks = [[] for i in range(0, len(freq_intv))]

# Loop over number_of_intervals
for i in range(0, number_of_intervals):
    # Pull current interval
    [p_temp, freq_temp] = [p_intv[i], freq_intv[i]]

    # Call peak finding function to get indices of peaks
    indx_peaks = detect_peaks(p_temp, mph=median_noise[i], show=False)

    # Save peak data in lists
    p_peaks[i] = [p_temp[x] for x in indx_peaks]
    freq_peaks[i] = [freq_temp[x] for x in indx_peaks]
    freq_muHz_peaks[i] = [freq_temp[x] / muHz for x in indx_peaks]


# # Calculate power/noise ratio
power_noise = as_f.func_over_noise(median_noise, p_peaks)


# # Finding peaks with power/noise > pn_limit

# Function call for power/noise check
[data_good, data_total] = as_f.good_peak_finder(power_noise, p_peaks, freq_peaks, freq_muHz_peaks, intv, pn_limit)
# Save peak data above power/noise limit (grouped so no intervals, peaks, freqs [muHz], freqs [Hz])
[p_g, f_m_g, f_g, pn_g] = data_good
# Save grouped total peak data, no intervals (total peaks, freq [muHz], freq [Hz])
[p_t, f_m_t, f_t] = data_total


# # Plot P/N and noise limit
plt.plot(f_m_t, power_noise, 'r*')
plt.plot(f_m_t, [pn_limit] * len(f_m_t), 'k')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.yscale('log', basey=10)
plt.ylabel('Power/noise ratio')
plt.legend(['Power Spectrum peaks', 'Chosen power/noise limit'])
plt.show()


# # Plot spectrum along with peaks over P/N limit
plt.plot(freq_muHz, p)
plt.plot(f_m_g, p_g, 'r*')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Power Spectrum')
plt.legend(['Power Spectrum', 'Peaks above power/noise limit'])
plt.show()


# # Find peaks that we believe in
repeat = input('Do you want to select intervals and peaks again (y/n)? If not, earlier selection will be used')
if repeat == 'y':
    plt.plot(freq_muHz, p)
    plt.plot(f_m_g, p_g, 'r*')
    plt.title('Select frequency intervals to view (Start+end for each)')
    plt.xlabel('Frequency [microHz]')
    plt.ylabel('Power Spectrum')
    plot_intervals = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    plt.close()

    intervals = []
    for i in range(0, len(plot_intervals)//2):
        [start, end] = [plot_intervals[2*i], plot_intervals[2*i+1]]
        [start, end] = [start[0], end[0]]
        plt.plot(freq_muHz, p)
        plt.plot(f_m_g, p_g, 'r*')
        plt.title('Select frequency intervals to view (Start+end for each)')
        plt.xlabel('Frequency [microHz]')
        plt.ylabel('Power Spectrum')
        plt.xlim([start, end])

        # Get y lim
        current_x_values_indx = np.where(np.logical_and(freq_muHz >= start, freq_muHz <= end))
        current_x_values_indx = list(current_x_values_indx)[0].tolist()
        current_y_values = [p[x] for x in current_x_values_indx]
        max_y_value = np.max(current_y_values)
        plt.ylim([-2, max_y_value+4])

        # Graphical input to select intervals
        intervals_temp = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
        plt.close()
        for k in range(0, len(intervals_temp)):
            intervals.append(intervals_temp[k])
            print(intervals)

    chosen_peak = []
    chosen_peak_freq = []
    for i in range(0, len(intervals)//2):
        [start, end] = [intervals[2*i], intervals[2*i+1]]
        [start, end] = [start[0], end[0]]
        current_interval = [start, end]

        # Find list indices for interval
        indices_p = np.where(np.logical_and(freq_muHz >= start, freq_muHz <= end))
        indices_peaks = np.where(np.logical_and(f_m_g >= start, f_m_g <= end))

        indices_p = list(indices_p)
        indices_p = indices_p[0].tolist()
        indices_peaks = list(indices_peaks)
        indices_peaks = indices_peaks[0].tolist()


        # Select values from lists

        freq_selected = [freq_muHz[x] for x in indices_p]
        p_selected = [p[x] for x in indices_p]
        p_g_selected = [p_g[x] for x in indices_peaks]
        f_g_selected = [f_m_g[x] for x in indices_peaks]

        # Plot interval
        plt.plot(freq_selected, p_selected)
        plt.plot(f_g_selected, p_g_selected, 'r*')
        plt.xlabel('Frequency [muHz]')
        plt.ylabel('Power')

        # Graphical input, select peaks you believe in
        plt.title('Select trustworthy peaks')
        trustworthy_peak_loc = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)

        # Find trustworthy peak from location
        chosen_peak_temp = []
        chosen_peak_freq_temp = []
        for k in range(0, len(trustworthy_peak_loc)):
            twpl_temp = trustworthy_peak_loc[k]
            twpl_temp = twpl_temp[0]
            twp_index = np.argmin([abs(twpl_temp - x) for x in f_g_selected])

            chosen_peak.append(p_g_selected[twp_index])
            chosen_peak_freq.append(f_g_selected[twp_index])
            chosen_peak_temp.append(p_g_selected[twp_index])
            chosen_peak_freq_temp.append(f_g_selected[twp_index])



        # Plot interval again to check
        plt.plot(freq_selected, p_selected)
        plt.plot(f_g_selected, p_g_selected, 'r*')
        plt.plot(chosen_peak_freq_temp, chosen_peak_temp, 'b*')
        plt.xlabel('Frequency [muHz]')
        plt.ylabel('Power')
        plt.show()

    ps_f.writer2col('found_peaks2', chosen_peak_freq, chosen_peak)

else:
    data = ps_f.dreader('found_peaks2.dat')
    chosen_peak_freq = []
    chosen_peak = []
    for row in data:
        chosen_peak_freq.append(row[0])
        chosen_peak.append(row[1])

plt.plot(freq_muHz, p)
plt.plot(f_m_g, p_g, 'r*')
plt.plot(chosen_peak_freq, chosen_peak, 'b*')
plt.xlabel('Frequency [muHz]')
plt.ylabel('Power')
#plt.yscale('log', basey=10)
plt.legend(['Power spectrum', 'Peaks found by peak finding', 'Trustworthy peaks'])
plt.show()

p_g_old = p_g
f_m_g_old = f_m_g

p_g = chosen_peak
f_m_g = chosen_peak_freq


if True:  # Fake block comment

    # # Autocorrelation with peaks only
    freq_mean_diff = np.mean(np.diff(freq_muHz))
    # Fill datasets with zeroes between peaks
    freq_filled = freq_muHz
    p_filled = [0] * len(freq_muHz)
    for i in range(0, len(p_g)):
        argument = [abs(f_m_g[i] - freq_muHz[x]) for x in range(0, len(freq_muHz))]
        indx = np.argmin(argument)
        # freq_filled[indx] = f_m_g[i]
        p_filled[indx] = p_g[i]

    # Use to full spectrum interval
    c1 = 0
    c2 = len(p_filled) - 1

    # Autocorrelation data select
    freq_ds2 = freq_filled[c1:c2]
    p_ds2 = p_filled[c1:c2]

    # Plot selected interval
    plt.plot(freq_filled, p_filled)
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Signal/Noise ratio')
    plt.legend(['S/N peak data'])
    plt.show()

    # Autocorrelation function call
    ac2 = as_f.autocorr(p_ds2)

    # Construct frequency steps

    freq_ac2 = np.linspace(0, freq_mean_diff * len(ac2), num=len(ac2), endpoint=False)

    # Limit to full frequency max
    cutindex = np.argmin([abs(x - freq_ac2[-1]) for x in freq_ac2])
    freq_ac2 = freq_ac2[0:cutindex]
    ac2 = ac2[0:cutindex]

    # Peak finding with specified min-peak-height
    ac2_peak_indx = detect_peaks(ac2, mph=0.02, show=False)
    ac2_peak = [ac2[x] for x in ac2_peak_indx]
    freq_ac2_peak = [freq_ac2[x] for x in ac2_peak_indx]

    # Plot autocorrelation
    plt.plot(freq_ac2, ac2)
    plt.plot(freq_ac2_peak, ac2_peak, 'r*')
    plt.xlabel('Frequency [microHz]')
    plt.xlim([-2, 400])
    plt.ylim([-0.01, 0.55])
    plt.ylabel('Autocorrelation')

    # Select autocorr frequencies
    indices = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    plt.close()
    large_freq = [None] * len(indices)
    ac2_large_freq = [None] * len(indices)
    for i in range(0, len(indices)):
        freq_chosen = indices[i]
        freq_chosen = freq_chosen[0]
        find_large_freq = np.argmin([abs(x - freq_chosen) for x in freq_ac2_peak])
        large_freq[i] = freq_ac2_peak[find_large_freq]
        ac2_large_freq[i] = ac2_peak[find_large_freq]

    large_frequency_split_list = 2 * np.diff(large_freq)
    large_frequency_split_mean = np.mean(large_frequency_split_list)
    large_frequency_split_std = np.std(large_frequency_split_list)
    print('Large frequency halfpoints')
    print(large_freq)
    print('Large frequency split', ', ', 'Standard deviation')
    print(large_frequency_split_mean, ', ', large_frequency_split_std)

    # Plot autocorrelation
    plt.plot(freq_ac2, ac2)
    plt.plot(freq_ac2_peak, ac2_peak, 'g.', markersize=7)
    plt.plot(large_freq, ac2_large_freq, 'r*', markersize=10)
    plt.legend(['Autocorrelation from power/noise peaks', 'Found autocorrelation peaks',
                'Selected large frequency separation peaks'])
    plt.xlim([-1, 400])
    plt.xlabel('Frequency [microHz]')
    plt.ylabel('Autocorrelation')
    plt.show()


# # Echelle diagram

# Input modulo value
modulus = 87.66

# Calculate logarithmic weights from peak height
base = 5
weights_peak = [math.log(x, base) for x in p_g]
weights_peak_old = [math.log(x, base) for x in p_g_old]
# weights_peak = [x ** 0.75 for x in p_g]
weights_full = [x ** 0.2 for x in p]
# weights_full = [math.log(x, base) for x in p]

# Call echelle2 with peak data to determine angular degree label placements
degree_coords = as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak,
                              xlabel='Frequency mod ' + str(modulus) + 'microHz',
                              ylabel='Frequency [microHz]', number_of_bins=100, mark_mode=True)

if True:
    # Call echelle2 with peak data, where no examination has occurred
    as_f.echelle2(f_m_g_old, modulus, heatmap=True, weights=weights_peak_old,
                  xlabel='Frequency mod ' + str(modulus) + 'microHz',
                  ylabel='Frequency [microHz]', number_of_bins=100, mark_mode=False, block=False)

    # Call echelle2 with peak data
    as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak,
                  xlabel='Frequency mod ' + str(modulus) + 'microHz',
                  ylabel='Frequency [microHz]', number_of_bins=100, mark_mode=False)

    # Call echelle2 with peak data to find small frequency separation
    points = as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak,
                           xlabel='Frequency mod ' + str(modulus) + 'microHz',
                           ylabel='Frequency [microHz]', number_of_bins=100, mark_mode=True,
                           degree_coords=degree_coords)

    # Get small frequency separation values from graphical input in echelle diagram
    small_freq_values = [None] * (len(points) // 2)
    frequencies_2 = [None] * (len(points) // 2)
    frequencies_0 = [None] * (len(points) // 2)
    mod_2 = [None] * (len(points) // 2)
    mod_0 = [None] * (len(points) // 2)
    for i in range(0, len(points) // 2):
        [point1, point2] = [points[2 * i], points[2 * i + 1]]
        x1, x2 = point1[0], point2[0]
        y1, y2 = point1[1], point2[1]
        small_freq_values[i] = x2 - x1
        frequencies_2[i], frequencies_0[i], mod_2[i], mod_0[i] = y1, y2, x1, x2

    small_freq_mean = np.mean(small_freq_values)
    small_freq_std = np.std(small_freq_values)

    # Print output
    print(' ')
    print('Echelle peak data')
    print('small_freq_values')
    print(small_freq_values)
    print('Small freq mean', ', ', 'Standard deviation')
    print(small_freq_mean, ', ', small_freq_std)

    # Print all frequencies and modulo
    print(' ')
    print('Echelle peak frequencies (l=2)')
    print(frequencies_2)
    print(' ')
    print('Echelle peak frequencies (l=0)')
    print(frequencies_0)
    print(' ')
    print('Echelle peak modulo (l=2)')
    print(mod_2)
    print(' ')
    print('Echelle peak modulo (l=0')
    print(mod_0)

    # Plot echelle diagram again, but with found points
    as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak,
                  xlabel='Frequency mod ' + str(modulus) + 'microHz', ylabel='Frequency [microHz]', number_of_bins=100,
                  mark_mode=False, degree_coords=degree_coords, plot_points=points)
