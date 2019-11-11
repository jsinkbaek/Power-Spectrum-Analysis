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


# # Frequency Conversion
muHz = 0.000001  # Variable for converting to microHz


# # Frequency calculation
resolution = 0.001 * muHz  # 0.01 normal
halfwidth = 6000 * muHz
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
ps_f.writer('spectrum', freq, spectral_power, res_pos.real, res_pos.imag)


# # Perform CLEAN procedure

# Create CLEAN window (index found by examining spectrum plot)
# window = range(1665000, 5271000)  # Sun signal_golf res 0.001muHz
window = range(195000, 500000)  # nuindi_filtered res 0.001muHz

# Call ps_f.clean_procedure
p_freq, p_power, p_a, p_b = ps_f.clean_procedure(time, flux, 50, halfwidth, resolution, window, mph=12)


# # Autocorrelation within window
acorr = as_f.autocorr(spectral_power[window])
freq_ac = np.arange(0, resolution*len(acorr), resolution)
print(type(acorr))
plt.figure()
plt.plot(freq_ac, acorr)
plt.show(block=False)

# # Peak autocorrelation
try:
    print(len(p_freq))
    print(len(p_power))
except NameError:
    data = ps_f.reader('clean_peaks')
    p_freq = data[0].ravel()
    p_power = data[1].ravel()

acorr_p, freq_ac_p = as_f.autocorr(p_power, number_of_steps=len(freq), x=p_freq, x_tot=freq)
plt.figure()
plt.plot(freq_ac_p, acorr_p)
plt.show(block=True)



