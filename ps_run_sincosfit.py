import ps_f
import matplotlib.pyplot as plt
# import cupy as cu
import numpy as np
import time as tm
import create_pspectrum
import clean_procedure
# import as_f
plt.rcParams.update({'font.size': 25})  # Changes font size for all plots

# # Pull variables from data
inpt1 = 'star2_corr'
# inpt1 = 'signal_betelgeuse'

data = ps_f.reader(inpt1)
time = data[0].ravel()  # (/ 86400 for betelgeuse)  # data[0] is in seconds
flux = data[1].ravel()


# # Frequency Conversion
muFreqUnit = 0.000001  # (* 86400 for betelgeuse)  # Variable for converting to microHz (or micro rev per day)

print('size:', flux.size)
# # Frequency calculation
resolution = 0.01 * muFreqUnit  # 0.01 normal
halfwidth = 500 * muFreqUnit  # (betelgeuse 20)
steps = int((2 * halfwidth) / resolution)

# # Spectrum calculation from sine cosine least squares fitting (betelgeuse 20.02 * muFreqUnit)
# results = create_pspectrum.cuda(flux, time, [6005 * muFreqUnit], halfwidth, resolution, chunk_size=100,
#                                dtype=cu.double)[0]
# results = create_pspectrum.numpy(flux, time, [(halfwidth+5) * muFreqUnit], halfwidth, resolution, chunk_size=500,
#                                 dtype=np.double)[0]
results = create_pspectrum.nufftpy(flux, time, half_width=2*halfwidth, resolution=resolution)
freq, spectral_power = results[0], results[1]

plt.plot(freq / muFreqUnit, spectral_power)
# plt.ylim([-10, 4*1.34*10**11])
# plt.xlim([0.04, 1.7])
plt.ylabel('Power')
plt.xlabel('muHz')
plt.legend(['sine/cosine fitting'])
plt.show()

# Estimate best resolution
resolution = ps_f.optimal_resolution(time, resolution, sfreq=halfwidth, half_width=10*muFreqUnit,
                                     chunk_size=2000)

# # Perform CLEAN procedure # #
# window = range(1665000, 5271000)
window = range(200000, 990000)
# mph = 0.00002
mph = 4
print('test')
p_freq, p_power, p_a, p_b = clean_procedure.numpy(time, flux, n_iter=50, freq_centre=(halfwidth+5*muFreqUnit),
                                                  resolution=resolution, window=window, mph=mph,
                                                  fft_half_width=halfwidth, chunk_size=2000, np_half_width=5*muFreqUnit)
try:
    print(len(p_freq))
    print(len(p_power))
except NameError:
    print('no peaks')
