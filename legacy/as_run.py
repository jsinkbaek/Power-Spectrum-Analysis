import math
import ps_f
import matplotlib.pyplot as plt
import numpy as np
from detect_peaks import detect_peaks
from scipy.optimize import curve_fit
import as_f
import time

plt.rcParams.update({'font.size': 22})


# # Data load from file
print('Input data to load, exclude extension (.dat)')
inpt1 = input()
data = ps_f.dreader(inpt1 + '.dat')

freq = []
p = []

for row in data:
    freq.append(row[0])
    p.append(row[1])
muHz = 0.000001
freq_muHz = [x / muHz for x in freq]

freq_mean_diff = np.mean(np.diff(freq_muHz))

# # Finding of data medians (noise level) from number_of_intervals different intervals (normal is 3)
[median_noise, intervals, p_intervals, freq_intervals] = as_f.median_noise(freq, p, number_of_intervals=3)

p_peaks = [[] for i in range(0, len(p_intervals))]
freq_peaks = [[] for i in range(0, len(freq_intervals))]
freq_muHz_peaks = [[] for i in range(0, len(freq_intervals))]


# # General peak finding (requirement that peak is larger than median), and calculation of peak/noise ratio
for i in range(0, len(p_intervals)):
    [p_intv, freq_intv] = [p_intervals[i], freq_intervals[i]]

    indices_peaks = detect_peaks(p_intv, mph=median_noise[i], show=False)
    p_peaks[i] = [p_intv[x] for x in indices_peaks]
    freq_peaks[i] = [freq_intv[x] for x in indices_peaks]
    freq_muHz_peaks[i] = [freq_intv[x] / muHz for x in indices_peaks]

peaks_over_noise = as_f.func_over_noise(median_noise, p_peaks)


# # Finding of peaks with power/noise > 6, based on the previous interval medians
[goodstuff, totalstuff] = as_f.good_peak_finder(peaks_over_noise, p_peaks, freq_peaks, freq_muHz_peaks, intervals)
[good_peaks, good_freq_muHz, good_freq] = goodstuff
[p_peaks_tot, freq_muHz_peaks_tot, freq_peaks_tot] = totalstuff


# # Plotting of spectrum along with "nice" peaks, where power/noise > 4
plt.plot(freq_muHz_peaks_tot, peaks_over_noise, 'r*')
plt.plot(freq_muHz_peaks_tot, [6] * len(freq_muHz_peaks_tot), 'k')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Power/noise ratio')
plt.legend(['Power Spectrum peaks', 'Chosen power/noise limit'])
plt.show()

#plt.plot(freq_muHz_peaks_tot, p_peaks_tot, 'r*')
#for i in range(0, len(intervals)):
#    plt.plot(freq_muHz_peaks[i], [4 * median_noise[i]] * len(freq_muHz_peaks[i]), 'k')

#plt.show()

plt.plot(freq_muHz, p)
plt.plot(good_freq_muHz, good_peaks, 'r*')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Power Spectrum')
plt.legend(['Power Spectrum', 'Peaks above power/noise limit'])
plt.show()


# # Recreating spectrum from before, just only with peaks and zeroes
freq_autocorr = freq_muHz
p_autocorr = [0] * len(freq_muHz)
for i in range(0, len(good_freq_muHz)):
    argument = [abs(good_freq_muHz[i] - freq_muHz[x]) for x in range(0, len(freq_muHz))]
    indx = np.argmin(argument)
    freq_autocorr[indx] = good_freq_muHz[i]
    p_autocorr[indx] = good_peaks[i]

# # Autocorrelation
plt.plot(p_autocorr)
plt.title('Select auto-correlation interval')
inpt2 = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
coord1 = inpt2[0]
coord1 = int(coord1[0])
coord2 = inpt2[1]
coord2 = int(coord2[0])
plt.close()

p_corr_select = p_autocorr[coord1:coord2]
freq_corr_select = freq_autocorr[coord1:coord2]

plt.plot(freq_autocorr, p_autocorr)
plt.plot(freq_corr_select, p_corr_select, 'r')
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Power Spectrum')
plt.legend(['Power spectrum reconstructed from peaks', 'Interval selected for autocorrelation'])
plt.show()

#  as_f.autocorrellator(freq_corr_select, p_corr_select)


correlation = as_f.autocorr(p_corr_select)

# # Reconstruct frequency steps
freq_autocorr_plot = np.linspace(0, freq_mean_diff*len(correlation), num=len(correlation), endpoint=False)
plt.plot(freq_autocorr_plot, correlation)
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Autocorrelation')

plt.show()

# # Find correlation frequencies
plt.plot(freq_autocorr_plot, correlation)
plt.xlim(-1, 0.7*freq_autocorr_plot[-1])
plt.ylim(-0.01, 0.4)
plt.xlabel('Frequency [' + r'mu' + 'Hz]')
plt.ylabel('Autocorrelation')
bigv = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
plt.close()
for i in range(0, len(bigv)):
    bigv_list = bigv[i]
    bigv[i] = bigv_list[0]
print(bigv)
deltav = 2 * np.mean(np.diff(bigv))
print(deltav)

# # Find index of our lower limit frequency (140 microHz)
indx_low = np.argmin([abs(x-140) for x in good_freq_muHz])
freq_echelle = good_freq_muHz[indx_low:]
p_echelle = good_peaks[indx_low:]
plt.plot(freq_echelle, p_echelle)
plt.show()

as_f.echelle(p_echelle, freq_echelle, deltav, block=True, weight=False)
#as_f.echelle(good_peaks, good_freq_muHz, deltav, block=True, weight=False)

deltav_found = 25
deltav = 20

startindx = np.argmin([abs(x-400) for x in freq_muHz])
endindx = np.argmin([abs(x-800) for x in freq_muHz])
while deltav < 1.5*deltav_found:
    as_f.echelle(p[startindx:endindx], freq_muHz[startindx:endindx], deltav, block=False, weight=True)
    deltav = deltav + 0.01
