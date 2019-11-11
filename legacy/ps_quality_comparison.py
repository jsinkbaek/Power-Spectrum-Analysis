"""
Script used to compare the median noise values for two different power spectra, by separating into similar noise
interval sections. It could f.ex. be used to compare power spectra obtained from the same data set subjected to
different filtering methods. Only used for an appendix section in the project.
"""

import matplotlib.pyplot as plt
import ps_f
import numpy as np
plt.rcParams.update({'font.size': 30})

n = 2
inpt = ['' for i in range(0, n)]
for i in range(0, n):
    print('Spectrum ' + str(i) + ', exclude extension')
    inpt[i] = input()


freq = [[] for i in range(0, n)]
freq_muHz = [[] for i in range(0, n)]
p = [[] for i in range(0, n)]

muHz = 0.000001

for i in range(0, n):
    data = ps_f.preader(inpt[i] + '.dat')
    for row in data:
        freq[i].append(row[0])
        freq_muHz[i].append(row[0] / muHz)
        p[i].append(row[1])


# # Interval selection
intervals = [[] for i in range(0, n)]
freq_intervals = [[] for i in range(0, n)]
h = 3

for i in range(0, n):
    plt.plot(p[i])
    plt.title('Select end/startpoint for first/second interval, end/startpoint for second/third interval')
    coords = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    coord1 = 0
    [coord2, coord3, coord4, coord5] = [coords[0], coords[0], coords[1], coords[1]]
    coord2 = int(coord2[0])
    coord3 = int(coord3[0])
    coord4 = int(coord4[0])
    coord5 = int(coord5[0])
    coord6 = len(p[i])
    intervals[i] = [[coord1, coord2], [coord3, coord4], [coord5, coord6]]
    if i == 1:
        interval_select = intervals[0]
        intervals[1] = intervals[0]
        [coord1, coord2] = interval_select[0]
        [coord3, coord4] = interval_select[1]
        [coord5, coord6] = interval_select[2]
    freq_temp = freq_muHz[i]
    freq_intervals[i] = [[freq_temp[coord1], freq_temp[coord2]], [freq_temp[coord3], freq_temp[coord4]],
                         [freq_temp[coord5], freq_temp[coord6-1]]]
    plt.close()

print('Freq Intervals:')
print(freq_intervals)

# #  Median of whole interval
median = [[] for i in range(0, n)]
colors = ['b', 'y', 'r']
a_colors = ['orange', 'purple', 'green']
for i in range(0, n):
    p_temp = p[i]
    intervals_temp = intervals[i]
    p_intervals = [[] for j in range(0, h)]
    median_intervals = [[] for j in range(0, h)]
    for j in range(0, h):
        interval_select = intervals_temp[j]
        p_intervals[j] = p_temp[interval_select[0]:interval_select[1]]
        median_intervals[j] = np.median(p_intervals[j])
    median[i] = [median_intervals[j] for j in range(0, len(median_intervals))]
    color = colors[i]
    print(inpt[i])
    print(median[i])

    freq_temp = freq_muHz[i]

    mead = median[i]
    legend = [''] * (2*h)
    for j in range(0, h):
        #plt.plot(range(0, 10), [mead[j]] * 10, color)
        interval_select = intervals_temp[j]
        freq_select = freq_temp[interval_select[0]:interval_select[1]]
        plt.plot(freq_select, p_intervals[j], colors[j])
        plt.plot(freq_select, [mead[j]] * len(freq_select), a_colors[j])
        legend[2*j] = 'Interval ' + str(j+1)
        legend[2*j+1] = 'Median ' + str(j+1)


    plt.xlabel('Frequency [microHz]')
    plt.ylabel('Power Spectrum')
    plt.legend(legend)
    plt.show()
