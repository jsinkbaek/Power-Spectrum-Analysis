import math
import ps_f
import matplotlib.pyplot as plt
import numpy as np
from detect_peaks import detect_peaks
from scipy.optimize import curve_fit
import as_f
import time
from pandas import Series

plt.rcParams.update({'font.size': 30})  # Changes font size for all plots

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

p_array = Series(p)
# # Smooth data
#p_smooth = as_f.gauss_smooth(p_array, window_size=99, sigma=20)
#p = p_smooth.tolist()

# # Find data medians (noise level) for all intervals
number_of_intervals = 3  # Sets the amount of medians necessary for the data (intervals)
[median_noise, intv, p_intv, freq_intv] = as_f.median_noise(freq, p, number_of_intervals=3)


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
pn_limit = 6

# Function call for power/noise check
[data_good, data_total] = as_f.good_peak_finder(power_noise, p_peaks, freq_peaks, freq_muHz_peaks, intv, pn_limit)
# Save peak data above power/noise limit (grouped so no intervals, peaks, freqs [muHz], freqs [Hz])
[p_g, f_m_g, f_g, pn_g] = data_good
# Save grouped total peak data, no intervals (total peaks, freq [muHz], freq [Hz])
[p_t, f_m_t, f_t] = data_total


# # Plot P/N and noise limit
plt.plot(f_m_t, power_noise, 'r*')
plt.plot(f_m_t, [pn_limit] * len(f_m_t), 'k')
#plt.plot(f_m_g, pn_g, 'b*')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.yscale('log', basey=10)
plt.ylabel('Power/noise ratio')
plt.legend(['Power Spectrum peaks', 'Chosen power/noise limit'])
plt.show()


# # Plot spectrum along with peaks over P/N limit
plt.plot(freq_muHz, p)
plt.plot(f_m_g, p_g, 'r*')
# 2plt.yscale('log', basey=10)
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Power Spectrum')
plt.legend(['Power Spectrum', 'Peaks above power/noise limit'])
plt.show()


# # Find expected oversampling peaks


# Plot and select peaks to examine
plt.plot(freq_muHz, p)
plt.plot(f_m_g, p_g, 'r*')
plt.title('Find oversampling peaks')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Power Spectrum')
plt.legend(['Power Spectrum', 'Peaks above power/noise limit'])
selected_frequencies = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
plt.close()

# Find the peaks selected
freq_index = [None] * len(selected_frequencies)
for i in range(0, len(selected_frequencies)):
    freq_temp = selected_frequencies[i]
    freq_temp = freq_temp[0]

    freq_index[i] = np.argmin([abs(freq_temp - x) for x in f_m_g])
freq_examine = [f_m_g[x] for x in freq_index]
peaks_examine = [p_g[x] for x in freq_index]


# Plot each with sinc function and mark oversampled peaks
limit = 20
f_m_g = np.asarray(f_m_g)
p_g = np.asarray(p_g)
for i in range(0, len(freq_examine)):
    indices = np.where(np.logical_and(f_m_g >= freq_examine[i] - limit, f_m_g <= freq_examine[i]+limit))

    sample_frequencies = [f_m_g[x] for x in indices]
    sample_peaks = [p_g[x] for x in indices]
    absolute_freq_difference = [abs(x-freq_examine[i]) for x in sample_frequencies]
    plt.figure()
    plt.plot(sample_frequencies, absolute_freq_difference, 'r*')
    plt.show(block=False)
    print(absolute_freq_difference)

    plt.figure()
    plt.plot(freq_muHz, p)
    plt.plot(f_m_g, p_g, 'r*')
    plt.xlim([freq_examine[i]-limit, freq_examine[i]+limit])
    plt.ylim([-2, peaks_examine[i]+2])
    plt.show(block=True)


if True:
    # # Autocorrelation with full dataset (in interval)

    ## Graphical input interval selection
    # plt.plot(p)
    # plt.title('Select autocorrelation interval')
    # inpt2 = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    # plt.close()
    ## Interval pull from input
    # [c1, c2] = [inpt2[0], inpt2[1]]
    # c1 = int(c1[0])  # x coordinate integer element
    # c2 = int(c2[0])

    # revert to full spectrum
    c1 = 0
    c2 = len(p) - 1

    # Autocorrelation data select
    freq_ds1 = freq_muHz[c1:c2]
    p_ds1 = p[c1:c2]

    # Plot selected interval
    plt.plot(freq_muHz, p)
    plt.plot(freq_ds1, p_ds1, 'r')
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Signal/Noise ratio')
    plt.legend(['S/N', 'Interval selected for autocorrelation'])
    plt.show()

    # Autocorrelation function call
    ac1 = as_f.autocorr(p_ds1)

    # Construct frequency steps
    freq_mean_diff = np.mean(np.diff(freq_muHz))
    freq_ac1 = np.linspace(0, freq_mean_diff * len(ac1), num=len(ac1), endpoint=False)

    # Peak finding
    ac1_peak_indx = detect_peaks(ac1, show=False)
    ac1_peak = [ac1[x] for x in ac1_peak_indx]
    freq_ac1_peak = [freq_ac1[x] for x in ac1_peak_indx]

    # Plot autocorrelation
    plt.plot(freq_ac1, ac1)
    plt.plot(freq_ac1_peak, ac1_peak, 'r*', markersize=7)
    plt.xlabel('Frequency [microHz]')
    plt.ylabel('Autocorrelation')
    plt.xlim([-1, freq_ac1[-1]])

    # Select autocorr frequencies for small frequency separation
    values = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    plt.close()

    small_freq_values = [None] * (len(values) // 2)
    plotpeaks_freq = [None] * len(values)
    plotpeaks_ac1 = [None] * len(values)
    for i in range(0, len(values) // 2):
        freq1, freq2 = values[2 * i], values[2 * i + 1]
        freq1, freq2 = freq1[0], freq2[0]
        find_nearest_peak_1 = np.argmin([abs(x - freq1) for x in freq_ac1_peak])
        find_nearest_peak_2 = np.argmin([abs(x - freq2) for x in freq_ac1_peak])

        small_freq_values[i] = freq_ac1_peak[find_nearest_peak_2] - freq_ac1_peak[find_nearest_peak_1]
        plotpeaks_freq[2 * i], plotpeaks_freq[2 * i + 1] = freq_ac1_peak[find_nearest_peak_1], freq_ac1_peak[
            find_nearest_peak_2]
        plotpeaks_ac1[2 * i], plotpeaks_ac1[2 * i + 1] = ac1_peak[find_nearest_peak_1], ac1_peak[find_nearest_peak_2]

    small_freq_mean = np.mean(small_freq_values)
    small_freq_std = np.std(small_freq_values)

    print('Small freq mean', ', ', 'Standard deviation')
    print(small_freq_mean, ', ', small_freq_std)

    # Plot autocorrelation
    plt.plot(freq_ac1, ac1)
    plt.plot(freq_ac1_peak, ac1_peak, 'g.', markersize=8)
    plt.plot(plotpeaks_freq, plotpeaks_ac1, 'r*', markersize=12)
    plt.legend(['Autocorrelation', 'Peaks', 'Peaks selected to find small frequency separation'])
    for i in range(0, len(plotpeaks_ac1) // 2):
        line = np.linspace(plotpeaks_freq[2 * i], plotpeaks_freq[2 * i + 1])
        y = [plotpeaks_ac1[2 * i]] * len(line)
        plt.plot(line, y, 'k', linewidth=0.8)

        endline1_y = np.linspace(plotpeaks_ac1[2 * i] - 0.01, plotpeaks_ac1[2 * i] + 0.01)
        endline1_x = [plotpeaks_freq[2 * i]] * len(endline1_y)
        plt.plot(endline1_x, endline1_y, 'k', linewidth=0.8)

        endline2_y = np.linspace(plotpeaks_ac1[2 * i] - 0.01, plotpeaks_ac1[2 * i] + 0.01)
        endline2_x = [plotpeaks_freq[2 * i + 1]] * len(endline2_y)
        plt.plot(endline2_x, endline2_y, 'k', linewidth=0.8)
    plt.legend(
        ['Autocorrelation', 'Peaks', 'Peaks selected to find small frequency separation', 'Small frequency separation'])
    plt.xlabel('Frequency [microHz]')
    plt.ylabel('Autocorrelation')
    plt.xlim([-1, freq_ac1[-1]])
    plt.show()

    # # Autocorrelation with peaks only

    # Fill datasets with zeroes between peaks
    freq_filled = freq_muHz
    p_filled = [0] * len(freq_muHz)
    for i in range(0, len(p_g)):
        argument = [abs(f_m_g[i] - freq_muHz[x]) for x in range(0, len(freq_muHz))]
        indx = np.argmin(argument)
        # freq_filled[indx] = f_m_g[i]
        p_filled[indx] = p_g[i]

    ## Graphical input interval selection
    # plt.plot(p_filled)
    # plt.title('Select autocorrelation interval')
    # inpt3 = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    # plt.close()

    ## Interval pull from input
    # [c1, c2] = [inpt3[0], inpt3[1]]
    # c1 = int(c1[0])  # x coordinate integer element
    # c2 = int(c2[0])

    # revert to full spectrum
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
    print('Large frequency split', ', ', 'Standard deviation')
    print(large_frequency_split_mean, ', ', large_frequency_split_std)

    # Plot autocorrelation
    plt.plot(freq_ac2, ac2)
    plt.plot(freq_ac2_peak, ac2_peak, 'g.', markersize=7)
    plt.plot(large_freq, ac2_large_freq, 'r*', markersize=10)
    plt.legend(['Autocorrelation from power/noise peaks', 'Found autocorrelation peaks',
                'Selected large frequency separation peaks'])
    plt.xlim([-1, freq_ac2[-1]])
    plt.xlabel('Frequency [microHz]')
    plt.ylabel('Autocorrelation')
    plt.show()


if False:
    # # Autocorrelation with full dataset (optional is specified interval)
    choose_interval = False

    if choose_interval is True:
        # Graphical input interval selection
        plt.plot(p)
        plt.title('Select autocorrelation interval')
        inpt2 = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
        plt.close()
        # Interval pull from input
        [c1, c2] = [inpt2[0], inpt2[1]]
        c1 = int(c1[0])  # x coordinate integer element
        c2 = int(c2[0])
    else:
        c1 = 0
        c2 = -1

    # Autocorrelation data select
    freq_ds1 = freq_muHz[c1:c2]
    p_ds1 = p[c1:c2]

    # Plot selected interval
    plt.plot(freq_muHz, p)
    plt.plot(freq_ds1, p_ds1, 'r')
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Signal/Noise ratio')
    plt.yscale('log', basey=10)
    plt.legend(['S/N', 'Interval selected for autocorrelation'])
    plt.show()

    # Autocorrelation function call
    ac1 = as_f.autocorr(p_ds1)

    # Construct frequency steps
    freq_mean_diff = np.mean(np.diff(freq_muHz))
    freq_ac1 = np.linspace(0, freq_mean_diff * len(ac1), num=len(ac1), endpoint=False)

    # Plot autocorrelation
    plt.plot(freq_ac1, ac1)
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Autocorrelation')
    plt.show()

    # # Autocorrelation with peaks only

    # Fill lists with zeroes between peaks

    # Preload lists
    f_fill = freq_muHz
    p_fill = [0] * len(freq_muHz)

    # Loop to fill peaks in p_fill
    for i in range(0, len(f_m_g)):
        argument = [abs(f_m_g[i] - freq_muHz[x]) for x in range(0, len(freq_muHz))]
        indx = np.argmin(argument)
        f_fill[indx] = f_m_g[i]  # Unnecessary
        p_fill[indx] = p_g[i]

    # Graphical input interval selection
    plt.plot(p_fill)
    plt.title('Select autocorrelation interval')
    inpt3 = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    plt.close()

    # Interval pull from input
    [c1, c2] = [inpt3[0], inpt3[1]]
    c1 = int(c1[0])  # x coordinate integer element
    c2 = int(c2[0])

    # Autocorrelation data select
    freq_ds2 = f_fill[c1:c2]
    p_ds2 = p_fill[c1:c2]

    # Plot selected interval
    plt.plot(f_fill, p_fill)
    plt.plot(freq_ds2, p_ds2, 'r')
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Signal/Noise ratio')
    plt.yscale('log', basey=10)
    plt.legend(['S/N', 'Interval selected for autocorrelation'])
    plt.show()

    # Autocorrelation function call
    ac2 = as_f.autocorr(p_ds2)

    # Construct frequency steps
    freq_ac2 = np.linspace(0, freq_mean_diff * len(ac2), num=len(ac2), endpoint=False)

    # Plot autocorrelation
    plt.plot(freq_ac2, ac2)
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Autocorrelation')
    plt.show()


# # Echelle diagram

# Input modulo value
modulus = input('What modulo value to use in Echelle diagram?')
modulus = float(modulus)

# Calculate logarithmic weights from peak height
base = 10
weights_peak = [math.log(x, base) for x in p_g]
# weights_peak = [x ** 0.75 for x in p_g]
# weights_full = [x ** 0.75 for x in p]
weights_full = [math.log(x, base) for x in p]

# Call echelle2 with peak data
degree_coords = as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak, xlabel='Frequency mod '+str(modulus)+'microHz',
              ylabel='Frequency [microHz]', number_of_bins=100, mark_mode=True)


# Call echelle2 with peak data to find small frequency separation
points = as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak, xlabel='Frequency mod '+str(modulus)+'microHz',
              ylabel='Frequency [microHz]', number_of_bins=100, ylim=[400, 1000], mark_mode=True, degree_coords=degree_coords)

# Get small frequency separation values from graphical input in echelle diagram
small_freq_values = [None] * (len(points)//2)
for i in range(0, len(points)//2):
    [point1, point2] = [points[2*i], points[2*i+1]]
    x1, x2 = point1[0], point2[0]
    small_freq_values[i] = x2 - x1
small_freq_mean = np.mean(small_freq_values)
small_freq_std = np.std(small_freq_values)

# Print output
print(' ')
print('Echelle peak data')
print('Small freq mean', ', ', 'Standard deviation')
print(small_freq_mean, ', ', small_freq_std)

# Plot echelle diagram again, but with found points
as_f.echelle2(f_m_g, modulus, heatmap=True, weights=weights_peak,
              xlabel='Frequency mod '+str(modulus)+'microHz',ylabel='Frequency [microHz]', number_of_bins=100,
              ylim=[400, 1000], mark_mode=False, degree_coords=degree_coords, plot_points=points)


# Call echelle2 with full data set
points = as_f.echelle2(freq_muHz, modulus, heatmap=True, weights=weights_full,
                       xlabel='Frequency mod '+str(modulus)+'microHz', ylabel='Frequency [microHz]',
                       number_of_bins=150, ylim=[400, 1000], degree_coords=degree_coords, mark_mode=True)

# Get small frequency separation values from graphical input in echelle diagram (full data set)
small_freq_values = [None] * (len(points)//2)
for i in range(0, len(points)//2):
    [point1, point2] = [points[2*i], points[2*i+1]]
    x1, x2 = point1[0], point2[0]
    small_freq_values[i] = x2 - x1
small_freq_mean = np.mean(small_freq_values)
small_freq_std = np.std(small_freq_values)

# Print output
print(' ')
print('Echelle all data')
print('Small freq mean', ', ', 'Standard deviation')
print(small_freq_mean, ', ', small_freq_std)

# Plot full echelle diagram with found points
as_f.echelle2(freq_muHz, modulus, heatmap=True, weights=weights_full, xlabel='Frequency mod '+str(modulus)+'microHz',
              ylabel='Frequency [microHz]', number_of_bins=150, ylim=[400, 1000], degree_coords=degree_coords,
              plot_points=points)

# Call echelle2 in interactive mode with peak data
repeat = True
stepsize = 0.5
print('')
mod = modulus
fig, ax = plt.subplots()
while repeat is True:
    plt.cla()
    [fig, ax] = as_f.echelle2(f_m_g, mod, heatmap=True, weights=weights_peak,
                              xlabel='Frequency mod ' + str(modulus) + 'microHz',
                              ylabel='Frequency [microHz]', animation_mode=True,
                              figure=[fig, ax], number_of_bins=100, ylim=[200, 450])
    plt.show(block=False)
    plt.pause(0.05)
    plt.close(2)
    key = input()
    if key == 'b':
        print('exiting')
        repeat = False
        break
    elif key == '+':
        stepsize = stepsize * 2
        print('Stepsize is now ' + str(stepsize))
    elif key == '-':
        stepsize = stepsize / 2
        print('Stepsize is now ' + str(stepsize))
    elif key == 'w':
        mod = mod + stepsize
        print('modulo is now ' + str(mod) + 'microHz')
    elif key == 's':
        mod = mod - stepsize
        print('modulo is now ' + str(mod) + 'microHz')


if False:
    # # Echelle diagram

    # Input modulo value
    modulus = input('What modulo value to use in Echelle diagram?')
    modulus = float(modulus)

    # Calculate logarithmic weights from peak height
    base = 2
    weights_peak = [math.log(x, base) for x in pn_g]
    # weights_peak = p_g
    # weights_peak = pn_g
    weights_full = [math.log(x, base) for x in p]

    # Call echelle2 with peak data
    as_f.echelle2(f_m_g, modulus, heatmap=False, weights=weights_peak,
                  xlabel='Frequency mod ' + str(modulus) + 'microHz',
                  ylabel='Frequency [microHz]', number_of_bins=60)

    # Call echelle2 with full data set
    as_f.echelle2(freq_muHz, modulus, heatmap=True, weights=weights_full,
                  xlabel='Frequency mod ' + str(modulus) + 'microHz',
                  ylabel='Frequency [microHz]', number_of_bins=60)

    # Call echelle2 in interactive mode with peak data
    repeat = True
    stepsize = 0.5
    print('')
    mod = modulus
    fig, ax = plt.subplots()
    while repeat is True:
        plt.cla()
        [fig, ax] = as_f.echelle2(f_m_g, mod, heatmap=True, weights=weights_peak,
                                  xlabel='Frequency mod ' + str(modulus) + 'microHz',
                                  ylabel='Frequency [microHz]', animation_mode=True,
                                  figure=[fig, ax], ylim=[400, 1000])
        plt.show(block=False)
        plt.pause(0.05)
        plt.close(2)
        key = input()
        if key == 'b':
            print('exiting')
            repeat = False
            break
        elif key == '+':
            stepsize = stepsize * 2
            print('Stepsize is now ' + str(stepsize))
        elif key == '-':
            stepsize = stepsize / 2
            print('Stepsize is now ' + str(stepsize))
        elif key == 'w':
            mod = mod + stepsize
            print('modulo is now ' + str(mod) + 'microHz')
        elif key == 's':
            mod = mod - stepsize
            print('modulo is now ' + str(mod) + 'microHz')

