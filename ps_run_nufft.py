import math
import ps_f
import matplotlib.pyplot as plt
import nufftpy
import numpy as np
import time as tm
import as_f

# # Input an call to load data set from file
print('What is the filename? (No file extension)')
inpt1 = input()

data = ps_f.reader(inpt1)

# # Pull variables from data
time = data[0].ravel()  # ravel() is used to change shape from (n, 1) to (n, )
flux = data[1].ravel()

# # Used to indicate which window, modulus and scaling should be used
dataset = 1  # 0 = sun, 1 = star2, 2 = nuindi

# # Frequency Conversion
muHz = 0.000001  # Variable for converting to microHz


# # Frequency calculation
resolution = 0.001 * muHz  # 0.01 normal
halfwidth_set = [6000, 1000, 600, 20]  # sun, star2, nuindi, betelgeuse
halfwidth = halfwidth_set[dataset] * muHz
steps = int((2 * halfwidth) / resolution)
freq = nufftpy.nufftfreqs(steps, df=resolution)
freq = freq[len(freq)//2:-1]


# # Spectrum calculation from non-uniform fft
result = nufftpy.nufft1(time, flux, steps, df=(resolution * 2 * math.pi))
res_pos = result[len(result)//2:-1]

spectral_power = res_pos.real ** 2 + res_pos.imag ** 2

plt.figure()
plt.plot(freq / muHz / 11.57, spectral_power)
# plt.plot(spectral_power)
plt.xlabel('Frequency [1/days]')
plt.xlim([0.04, 1.7])
plt.ylim([-10, 1.34*10**11])
plt.ylabel('Power (non-uniform fft)')
plt.show()

# ps_f.writer('spectrum', freq, spectral_power, res_pos.real, res_pos.imag)


# # Perform CLEAN procedure

# Create CLEAN window (index found by examining spectrum plot)
window_set = [range(1665000, 5271000), range(200000, 990000), range(195000, 500000), None]  # sun, star2, nuindi, betty
window = window_set[dataset]

# Call ps_f.clean_procedure
mph_set = [0.0002, 4, 12, 0.2 * 10**9]  # sun, star2, nuindi, betelgeuse
mph = mph_set[dataset]
p_freq, p_power, p_a, p_b = ps_f.clean_procedure(time, flux, 200, halfwidth, resolution, window, mph=mph)

try:
    print(len(p_freq))
    print(len(p_power))
except NameError:
    data = ps_f.reader('clean_peaks')
    p_freq = data[0].ravel()
    p_power = data[1].ravel()

# # Autocorrelation within window
if window is not None:
    acorr = as_f.autocorr(spectral_power[window])
else:
    acorr = as_f.autocorr(spectral_power)
freq_ac = np.arange(0, resolution * len(acorr), resolution)
plt.figure()
plt.plot(freq_ac, acorr)
plt.show(block=False)

# # Peak autocorrelation
acorr_p, freq_ac_p = as_f.autocorr(p_power, number_of_steps=len(freq), x=p_freq, x_tot=freq)
plt.figure()
plt.plot(freq_ac_p, acorr_p)
plt.show(block=True)


# # Echelle diagram
p_freq = np.asarray(p_freq)
weights = np.sqrt(p_power)
modulus_set = [136, 87.66, 25.21775, 3.43]  # sun, star2, nuindi, betelgeuse
modulus = modulus_set[dataset]
as_f.echelle(p_freq / muHz, modulus, heatmap=True, weights=weights, number_of_bins=64)


