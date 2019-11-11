import matplotlib.pyplot as plt
import ps_f
plt.rcParams.update({'font.size': 30})

print('How many sets do you want to plot? Up to 2.')
inpt0 = input()
print('Which spectrum to plot? Insert datafile, excluding .dat')
inpt1 = input()

freq = []
p = []
muHz = 0.000001

if inpt1 == 'corrected_nuindi120s':
    data = ps_f.dreader(inpt1 + '.dat')
else:
    data = ps_f.dreader(inpt1 + '.dat')
for row in data:
    freq.append(row[0])
    p.append(row[1])

if inpt0 == 'data_used_to_fit_md':
    freq_muHz = freq
else:
    freq_muHz = [x / muHz for x in freq]



plt.subplot(1, 1, 1)
plt.plot(freq_muHz, p, 'b.')
plt.plot(freq_muHz, p, 'b', linewidth=0.2)
#plt.yscale('log', basey=10)
#plt.xscale('log', basex=10)
#plt.xlim([0, 1000])
#plt.ylim([0, 1150])
plt.ylabel('Power Spectrum')
plt.xlabel('Frequency [microHz]')

if inpt0 == '2':
    print('Which second spectrum to plot?')
    inpt2 = input()
    freq1 = []
    p1 = []
    data1 = ps_f.dreader(inpt2 + '.dat')
    for row in data1:
        freq1.append(row[0])
        p1.append(row[1])
    if inpt1 == 'data_used_to_fit_nd':
        freq1_muHz = freq1
    else:
        freq1_muHz = [x / muHz for x in freq1]
    plt.subplot(2, 1, 2)
    plt.plot(freq1_muHz, p1)
    plt.xlabel('Frequency [microHz]')
    #plt.yscale('log', basey=10)
    #plt.xscale('log', basex=10)
    #plt.xlim([0, 1000])
    #plt.ylim([0, 1150])
    plt.ylabel('Power Spectrum')

plt.show()
