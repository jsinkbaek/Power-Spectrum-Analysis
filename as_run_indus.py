"""
Script to run through all (except P/N calculation, done in matlab) the power-spectrum analysis for nu Indi.
-   Does peak-finding by using the function detect_peaks (https://github.com/demotu/BMC). Assumes that the spectrum
    imported from dat-file is a P/N spectrum when doing this.
-   Calculates autocorrelation using the whole P/N data-set, in order to (using graphical input selection) determine
    small frequency separation.
-   Calculates autocorrelation using only data from peaks in P/N data-set, in order to (using graphical input selection)
    determine large frequency separation.
-   Plots echelle heatmap diagrams multiple times using only data from peaks in P/N data-set, in order to determine
    small frequency separation.
-   Plots echelle heatmap diagrams multiple times using the whole P/N data-set, in order to determine small frequency
    separation (again)
-   Plots echelle heatmap diagram in animation mode. Not used for anything specific.
"""

import ps_f
import matplotlib.pyplot as plt
import numpy as np
from detect_peaks import detect_peaks
import as_f

plt.rcParams.update({'font.size': 22})  # Changes font size for all plots

# # Data load from file
# Input data file name
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


# # Data plot and interval select
# plt.plot(p)
# [coord1, coord2] = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
# plt.close()
# [coord1, coord2] = [int(coord1[0]), int(coord2[0])]
[coord1, coord2] = [13777, 57688]  # Hardcoded for repetition, used commented out part to determine
freq_old = freq_muHz
freq_muHz = freq_muHz[coord1:coord2]

p_old = p
p = p[coord1:coord2]

# Plot selection
plt.plot(freq_old, p_old)
plt.plot(freq_muHz, p, 'r')
plt.xlabel('Frequency [microHz]')
plt.ylabel('Power/Noise')
plt.legend(['Power/Noise spectrum', 'Region selected for analysis'])
plt.show()

# # Peak finding
# Assumes that p is a signal-to-noise ratio already. Otherwise use as_run.py

# Function call, peak index finder, requires min. 6 signal/noise ratio
sn_limit = 6
indices_peaks = detect_peaks(p, mph=sn_limit, show=False)

# Peak values and frequencies found
p_peaks = [p[x] for x in indices_peaks]
freq_peaks = [freq[x] for x in indices_peaks]
freq_muHz_peaks = [freq_muHz[x] for x in indices_peaks]


# # Limit plotting and peaks
plt.plot(freq_muHz, p)
plt.plot(freq_muHz, [sn_limit] * len(freq_muHz), 'k')
plt.plot(freq_muHz_peaks, p_peaks, 'r*', markersize=8)
plt.xlabel('Frequency [microHz]')
plt.ylabel('Power/noise ratio')
plt.legend(['Power/noise Spectrum', 'Chosen power/noise limit', 'Peaks above power/noise limit'])
plt.show()


# # Peak and signal to noise ratio plotting
plt.plot(freq_muHz, p)
plt.plot(freq_muHz_peaks, p_peaks, 'r*')
plt.xlabel('Frequency [microHz]')
plt.ylabel('Power/Noise ratio')
plt.legend(['Power/Noise data', 'Peaks above power/noise limit'])
plt.show()


# # # # # # # # # # # # # # # Autocorrelation # # # # # # # # # # # # # # # #


# # # # Autocorrelation with full dataset (in interval) # # # #

# Use to full spectrum
c1 = 0
c2 = len(p)-1

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
freq_ac1 = np.linspace(0, freq_mean_diff*len(ac1), num=len(ac1), endpoint=False)

# Peak finding
ac1_peak_indx = detect_peaks(ac1, show=False)
ac1_peak = [ac1[x] for x in ac1_peak_indx]
freq_ac1_peak = [freq_ac1[x] for x in ac1_peak_indx]

# Plot autocorrelation
plt.plot(freq_ac1, ac1)
plt.plot(freq_ac1_peak, ac1_peak, 'r*', markersize=7)
plt.xlabel('Frequency [microHz]')
plt.ylabel('Autocorrelation')
plt.xlim([-1, 110])

# Select autocorr frequencies for small frequency separation
values = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
plt.close()

# Pre-load lists for small frequency separation and plot-points
small_freq_values = [None] * (len(values)//2)
plotpeaks_freq = [None] * len(values)
plotpeaks_ac1 = [None] * len(values)

# Loop through to fill lists with data from values variable
for i in range(0, len(values)//2):
    # Pull frequency data from values variable
    freq1, freq2 = values[2*i], values[2*i+1]
    freq1, freq2 = freq1[0], freq2[0]

    # Find the selected peak indices from graphical input frequencies
    find_nearest_peak_1 = np.argmin([abs(x - freq1) for x in freq_ac1_peak])
    find_nearest_peak_2 = np.argmin([abs(x - freq2) for x in freq_ac1_peak])

    # Calculate small frequency separation from selected peak frequencies
    small_freq_values[i] = freq_ac1_peak[find_nearest_peak_2] - freq_ac1_peak[find_nearest_peak_1]

    # Save selected peak frequency and autocorrelation points for plotting
    plotpeaks_freq[2*i], plotpeaks_freq[2*i+1] = freq_ac1_peak[find_nearest_peak_1], freq_ac1_peak[find_nearest_peak_2]
    plotpeaks_ac1[2*i], plotpeaks_ac1[2*i+1] = ac1_peak[find_nearest_peak_1], ac1_peak[find_nearest_peak_2]

# Calculate small frequency separation mean and standard deviation
small_freq_mean = np.mean(small_freq_values)
small_freq_std = np.std(small_freq_values)

# Print
print('Small freq values')
print(small_freq_values)
print('Small freq mean', ', ', 'Standard deviation')
print(small_freq_mean, ', ', small_freq_std)

# Plot autocorrelation with peaks selected
plt.plot(freq_ac1, ac1)
plt.plot(freq_ac1_peak, ac1_peak, 'g.', markersize=8)
plt.plot(plotpeaks_freq, plotpeaks_ac1, 'r*', markersize=12)
plt.legend(['Autocorrelation', 'Peaks', 'Peaks selected to find small frequency separation'])
for i in range(0, len(plotpeaks_ac1)//2):
    line = np.linspace(plotpeaks_freq[2*i], plotpeaks_freq[2*i+1])
    y = [plotpeaks_ac1[2*i]] * len(line)
    plt.plot(line, y, 'k', linewidth=0.8)

    endline1_y = np.linspace(plotpeaks_ac1[2*i]-0.01, plotpeaks_ac1[2*i]+0.01)
    endline1_x = [plotpeaks_freq[2*i]] * len(endline1_y)
    plt.plot(endline1_x, endline1_y, 'k', linewidth=0.8)

    endline2_y = np.linspace(plotpeaks_ac1[2 * i] - 0.01, plotpeaks_ac1[2 * i] + 0.01)
    endline2_x = [plotpeaks_freq[2 * i + 1]] * len(endline2_y)
    plt.plot(endline2_x, endline2_y, 'k', linewidth=0.8)
plt.legend(['Autocorrelation', 'Peaks', 'Peaks selected to find small frequency separation',
            'Small frequency separation'])
plt.xlabel('Frequency [microHz]')
plt.ylabel('Autocorrelation')
plt.xlim([-1, 110])
plt.show()


# # # # Autocorrelation with peaks only # # # #

# Fill datasets with zeroes between peaks
freq_filled = freq_muHz
p_filled = [0] * len(freq_muHz)
for i in range(0, len(p_peaks)):
    argument = [abs(freq_muHz_peaks[i]-freq_muHz[x]) for x in range(0, len(freq_muHz))]
    indx = np.argmin(argument)
    p_filled[indx] = p_peaks[i]

# Use whole interval
c1 = 0
c2 = len(p_filled)-1

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
freq_ac2 = np.linspace(0, freq_mean_diff*len(ac2), num=len(ac2), endpoint=False)

# Limit to 110 frequency max
cutindex = np.argmin([abs(x-110) for x in freq_ac2])
freq_ac2 = freq_ac2[0:cutindex]
ac2 = ac2[0:cutindex]

# Peak finding with specified min-peak-height
ac2_peak_indx = detect_peaks(ac2, mph=0.04, show=False)
ac2_peak = [ac2[x] for x in ac2_peak_indx]
freq_ac2_peak = [freq_ac2[x] for x in ac2_peak_indx]

# Plot autocorrelation
plt.plot(freq_ac2, ac2)
plt.plot(freq_ac2_peak, ac2_peak, 'r*')
plt.xlabel('Frequency [microHz]')
plt.ylabel('Autocorrelation')


# # Select autocorr frequencies
indices = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
plt.close()

# Preload lists for selected peaks data
large_freq = [None] * len(indices)
ac2_large_freq = [None] * len(indices)

# Loop to find selected peaks from graphical input coordinates
for i in range(0, len(indices)):
    freq_chosen = indices[i]
    freq_chosen = freq_chosen[0]
    find_large_freq = np.argmin([abs(x - freq_chosen) for x in freq_ac2_peak])
    large_freq[i] = freq_ac2_peak[find_large_freq]
    ac2_large_freq[i] = ac2_peak[find_large_freq]

# Determine large frequency separation as difference in frequency between selected peaks
large_frequency_split_list = 2 * np.diff(large_freq)
large_frequency_split_mean = np.mean(large_frequency_split_list)
large_frequency_split_std = np.std(large_frequency_split_list)
print('Large frequency half-values')
print(large_freq)
print('Large frequency split', ', ', 'Standard deviation')
print(large_frequency_split_mean, ', ', large_frequency_split_std)

# Plot autocorrelation
plt.plot(freq_ac2, ac2)
plt.plot(freq_ac2_peak, ac2_peak, 'g.', markersize=7)
plt.plot(large_freq, ac2_large_freq, 'r*', markersize=10)
plt.legend(['Autocorrelation from power/noise peaks', 'Found autocorrelation peaks',
            'Selected large frequency separation peaks'])
plt.xlim([-1, 110])
plt.xlabel('Frequency [microHz]')
plt.ylabel('Autocorrelation')
plt.show()


# # # # # # # # # # # # Echelle diagram # # # # # # # # # # # # # # # #

# Input modulo value
modulus = input('What modulo value to use in Echelle diagram?')
modulus = float(modulus)

# Calculate logarithmic weights from peak height
base = 2
# weights_peak = [math.log(x, base) for x in p_peaks]
weights_peak = [x ** 0.75 for x in p_peaks]
weights_full = [x ** 0.75 for x in p]
# weights_full = [math.log(x, base) for x in p]

# Call echelle2 with peak data
degree_coords = as_f.echelle2(freq_muHz_peaks, modulus, heatmap=True, weights=weights_peak,
                              xlabel='Frequency mod '+str(modulus)+'microHz', ylabel='Frequency [microHz]',
                              number_of_bins=100, ylim=[250, 400], mark_mode=True)


# Call echelle2 with peak data to find small frequency separation
points = as_f.echelle2(freq_muHz_peaks, modulus, heatmap=True, weights=weights_peak,
                       xlabel='Frequency mod '+str(modulus)+'microHz',
                       ylabel='Frequency [microHz]', number_of_bins=100, ylim=[200, 450],
                       mark_mode=True, degree_coords=degree_coords)

# Get small frequency separation values from graphical input in echelle diagram
small_freq_values = [None] * (len(points)//2)
frequencies_2 = [None] * (len(points) // 2)
frequencies_0 = [None] * (len(points) // 2)
mod_2 = [None] * (len(points) // 2)
mod_0 = [None] * (len(points) // 2)
for i in range(0, len(points)//2):
    [point1, point2] = [points[2*i], points[2*i+1]]
    x1, x2 = point1[0], point2[0]
    y1, y2 = point1[1], point2[1]
    small_freq_values[i] = x2 - x1
    frequencies_2[i], frequencies_0[i], mod_2[i], mod_0[i] = y1, y2, x1, x2

small_freq_mean = np.mean(small_freq_values)
small_freq_std = np.std(small_freq_values)

# Print output
print(' ')
print('Echelle peak data')
print('Small frequency values')
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
as_f.echelle2(freq_muHz_peaks, modulus, heatmap=True, weights=weights_peak,
              xlabel='Frequency mod '+str(modulus)+'microHz',ylabel='Frequency [microHz]', number_of_bins=100,
              ylim=[200, 450], mark_mode=False, degree_coords=degree_coords, plot_points=points)


# Call echelle2 with full data set
points = as_f.echelle2(freq_muHz, modulus, heatmap=True, weights=weights_full,
                       xlabel='Frequency mod '+str(modulus)+'microHz', ylabel='Frequency [microHz]',
                       number_of_bins=150, ylim=[200, 450], degree_coords=degree_coords, mark_mode=True)

# Get small frequency separation values from graphical input in echelle diagram (full data set)
small_freq_values = [None] * (len(points)//2)
frequencies_2 = [None] * (len(points) // 2)
frequencies_0 = [None] * (len(points) // 2)
mod_2 = [None] * (len(points) // 2)
mod_0 = [None] * (len(points) // 2)
for i in range(0, len(points)//2):
    [point1, point2] = [points[2*i], points[2*i+1]]
    x1, x2 = point1[0], point2[0]
    y1, y2 = point1[1], point2[1]
    small_freq_values[i] = x2 - x1
    frequencies_2[i], frequencies_0[i], mod_2[i], mod_0[i] = y1, y2, x1, x2

small_freq_mean = np.mean(small_freq_values)
small_freq_std = np.std(small_freq_values)

# Print output
print(' ')
print('Echelle all data')
print('small_freq_values')
print(small_freq_values)
print('Small freq mean', ', ', 'Standard deviation')
print(small_freq_mean, ', ', small_freq_std)

# Print all frequencies and modulo
print(' ')
print('Echelle frequencies (l=2)')
print(frequencies_2)
print(' ')
print('Echelle frequencies (l=0)')
print(frequencies_0)
print(' ')
print('Echelle modulo (l=2)')
print(mod_2)
print(' ')
print('Echelle modulo (l=0')
print(mod_0)

# Plot full echelle diagram with found points
as_f.echelle2(freq_muHz, modulus, heatmap=True, weights=weights_full, xlabel='Frequency mod '+str(modulus)+'microHz',
              ylabel='Frequency [microHz]', number_of_bins=150, ylim=[200, 450], degree_coords=degree_coords,
              plot_points=points)

# Call echelle2 in interactive mode with peak data
repeat = True
stepsize = 0.5
print('')
mod = modulus
fig, ax = plt.subplots()
while repeat is True:
    plt.cla()
    [fig, ax] = as_f.echelle2(freq_muHz_peaks, mod, heatmap=True, weights=weights_peak,
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


