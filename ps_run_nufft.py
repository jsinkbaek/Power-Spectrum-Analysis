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
dataset = 1  # 1 = sun, 2 = star2, 3 = nuindi

# # Frequency Conversion
muHz = 0.000001  # Variable for converting to microHz


# # Frequency calculation
resolution = 0.001 * muHz  # 0.01 normal
halfwidth_set = [6000, 1000, 600]  # sun, star2, nuindi
# halfwidth = 6000 * muHz  # sun
halfwidth = halfwidth_set[dataset] * muHz
steps = int((2 * halfwidth) / resolution)
freq = nufftpy.nufftfreqs(steps, df=resolution)
freq = freq[len(freq)//2:-1]


# # Spectrum calculation from non-uniform fft
result = nufftpy.nufft1(time, flux, steps, df=(resolution * 2 * math.pi))
res_pos = result[len(result)//2:-1]

spectral_power = res_pos.real ** 2 + res_pos.imag ** 2
plt.plot(spectral_power)
# plt.plot(res_pos.real, 'r--', linewidth=0.3)
# plt.plot(res_pos.imag, 'b--', linewidth=0.3)
plt.show()
# ps_f.writer('spectrum', freq, spectral_power, res_pos.real, res_pos.imag)


# # Perform CLEAN procedure

# Create CLEAN window (index found by examining spectrum plot)
window_set = [range(1665000, 5271000), range(8400, 990000), range(195000, 500000)]  # sun, star2, nuindi
#window = range(1665000, 5271000)  # Sun signal_golf res 0.001muHz
window = window_set[dataset]

# Call ps_f.clean_procedure
mph_set = [0.0002, 2, 12]  # sun, star2, nuindi
# mph = 0.0002  # golf (sun)
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
acorr = as_f.autocorr(spectral_power[window])
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
weights = np.log2(p_power)
modulus_set = [136, 87.66, 25.21775]  # sun, star2, nuindi
# modulus = 136  # sun
modulus = modulus_set[dataset]
as_f.echelle(p_freq / muHz, modulus, heatmap=True, weights=weights, number_of_bins=124)


