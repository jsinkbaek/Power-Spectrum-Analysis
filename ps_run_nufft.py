import math
import ps_f
import matplotlib.pyplot as plt
import nufftpy
import numpy as np
import time as tm

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
t1 = tm.time()
resolution = 0.001 * muHz  # 0.01 normal
halfwidth = 6000 * muHz
steps = int((2 * halfwidth) / resolution)
freq = nufftpy.nufftfreqs(steps, df=resolution)
freq = freq[len(freq)//2:-1]
t2 = tm.time()
print(t2-t1)

# # Spectrum calculation from Non-uniform fft
result = nufftpy.nufft1(time, flux, steps, df=(resolution * 2 * math.pi))
t3 = tm.time()
print(t3-t2)
res_pos = result[len(result)//2:-1]

spectral_power = res_pos.real ** 2 + res_pos.imag ** 2

plt.plot(spectral_power)
# plt.plot(res_pos.real, 'r--', linewidth=0.3)
# plt.plot(res_pos.imag, 'b--', linewidth=0.3)
plt.show()
t4 = tm.time()
print(t4-t3)
# inpt2 = input('Name of file to save spectrum in?')
ps_f.writer('spectrum', freq, spectral_power, res_pos.real, res_pos.imag)

window = range(1665000, 5271000)

ps_f.clean_procedure(time, flux, 400, halfwidth, resolution, window, mph=0.00002)
