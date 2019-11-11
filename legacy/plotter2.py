import matplotlib.pyplot as plt
import ps_f
import numpy as np
plt.rcParams.update({'font.size': 25})


inpt1 = 'data_used_to_fit_md'
inpt01 = 'md_ps'
inpt11 = 'md_fit'

inpt2 = 'md_sn'



freq = []
p = []
muHz = 0.000001

data = ps_f.dreader(inpt1 + '.dat')

for row in data:
    freq.append(row[0])
    p.append(row[1])


low_index = np.argmin([abs(x-226) for x in freq])
high_index = np.argmin([abs(x-375) for x in freq])

freq_standard = []
p_standard = []
data = ps_f.dreader(inpt01 + '.dat')
for row in data:
    freq_standard.append(row[0])
    p_standard.append(row[1])
freq_standard[:] = [x / muHz for x in freq_standard]

data01 = ps_f.dreader(inpt11 + '.dat')
fit = []
freq01 = []
for row in data01:
    freq01.append(row[0])
    fit.append(row[1])

freq01_muHz = [x / muHz for x in freq01]
plt.subplot(2, 2, 1)
plt.plot(freq_standard, p_standard, 'C1--', linewidth=1)
plt.plot(freq[0:low_index+1], p[0:low_index+1], 'C0')
plt.plot(freq[high_index:len(freq)+1], p[high_index:len(p)+1], 'C0')
plt.plot(freq01_muHz, fit, 'r')
plt.yscale('log', basey=10)
plt.xscale('log', basex=10)
#plt.xlim([0, 1000])
plt.ylim([0.0001, 1500])
plt.ylabel('Power Spectrum')

freq1 = []
p1 = []
data1 = ps_f.dreader(inpt2 + '.dat')
for row in data1:
    freq1.append(row[0])
    p1.append(row[1])

freq1[:] = [x / muHz for x in freq1]

plt.subplot(2, 2, 3)
plt.plot(freq1, p1)
plt.xlabel('Frequency [microHz]')
plt.yscale('log', basey=10)
plt.xscale('log', basex=10)
# plt.xlim([0, 1000])
plt.ylim([0.0001, 1500])
plt.ylabel('Signal/Noise')


plt.subplot(2, 2, 2)
plt.plot(freq_standard, p_standard, 'C1--', linewidth=1)
plt.plot(freq[0:low_index+1], p[0:low_index+1], 'C0')
plt.plot(freq[high_index:len(freq)+1], p[high_index:len(p)+1], 'C0')
plt.plot(freq01_muHz, fit, 'r')
#plt.yscale('log', basey=10)
#plt.xscale('log', basex=10)
plt.xlim([0, 1000])
plt.ylim([-1, 1200])
#plt.ylabel('Power Spectrum')


plt.subplot(2, 2, 4)
plt.plot(freq1, p1)
plt.xlabel('Frequency [microHz]')
#plt.yscale('log', basey=10)
#plt.xscale('log', basex=10)
plt.xlim([0, 1000])
#plt.ylim([-1, 1200])
#plt.ylabel('Signal/Noise')

plt.show()