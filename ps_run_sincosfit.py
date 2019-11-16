import math
import ps_f
import matplotlib.pyplot as plt
import numpy as np
import time as tm
# import as_f

# # Pull variables from data
inpt1 = 'signal_golf'

data = ps_f.reader(inpt1)
time = data[0].ravel()
flux = data[1].ravel()

# # Frequency Conversion
muHz = 0.000001  # Variable for converting to microHz

# # Frequency calculation
resolution = 0.01 * muHz  # 0.01 normal
halfwidth = 6000*muHz
steps = int((2 * halfwidth) / resolution)

# # Spectrum calculation from sine cosine least squares fitting
results = ps_f.create_pspectrum(flux, time, [6005*muHz], halfwidth, resolution, chunk_size=800, dtype=np.longdouble)[0]
freq, spectral_power = results[0], results[1]

plt.plot(freq, spectral_power)
plt.show()
