import math
import ps_f
import matplotlib.pyplot as plt
import nufftpy
import numpy as np

# # Input an call to load data set from file
print('What is the filename? (No file extension)')
inpt1 = input()

data = ps_f.reader(inpt1 + '.dat')

# # Preload variables to append data from file to
tid = []
cflux = []
mcflux = []
variance = []

# # Loop to append data to variables
for row in data:
    tid.append(row[0])
    cflux.append(row[1])
    mcflux.append(row[2])
    variance.append(row[3])

# # Frequency Conversion
muHz = 0.000001  # Variable for converting to microHz

# # Frequency calculation
resolution = 0.01 * muHz  # 0.01 normal
halfwidth = 6000 * muHz
steps = int((2 * halfwidth) / resolution)
freq = nufftpy.nufftfreqs(steps, df=resolution)
freq = freq[len(freq)//2:-1]

# # Spectrum calculation from Non-uniform fft

result = nufftpy.nufft1(tid, cflux, steps, df=(resolution * 2 * math.pi))
result = result[len(result)//2:-1]

spectral_power = result.real ** 2 + result.imag ** 2

plt.plot(freq, spectral_power)
plt.show()

inpt2 = input('Name of file to save spectrum in?')
ps_f.pwriter(inpt2, freq, spectral_power, result.real, result.imag)

