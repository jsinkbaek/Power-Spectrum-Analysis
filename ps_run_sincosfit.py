import ps_f
import matplotlib.pyplot as plt
import numpy as np
import time as tm
# import as_f
plt.rcParams.update({'font.size': 25})  # Changes font size for all plots

# # Pull variables from data
# inpt1 = 'signal_golf'
inpt1 = 'signal_betelgeuse'

data = ps_f.reader(inpt1)
time = data[0].ravel() / 86400  # data[0] is in seconds
flux = data[1].ravel()

# # Frequency Conversion
muFreqUnit = 0.000001 * 86400  # Variable for converting to microHz (or micro rev per day)

# # Frequency calculation
resolution = 0.01 * muFreqUnit  # 0.01 normal
halfwidth = 20 * muFreqUnit
steps = int((2 * halfwidth) / resolution)

# # Spectrum calculation from sine cosine least squares fitting
results = ps_f.create_pspectrum(flux, time, [20.02 * muFreqUnit], halfwidth, resolution, chunk_size=50,
                                dtype=np.longdouble)[0]
freq, spectral_power = results[0], results[1]

plt.plot(freq, spectral_power)
plt.ylim([-10, 4*1.34*10**11])
plt.xlim([0.04, 1.7])
plt.ylabel('Power')
plt.xlabel('1/days')
plt.legend(['sine/cosine fitting'])
plt.show()
