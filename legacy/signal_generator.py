from __future__ import division, print_function

import matplotlib.pyplot as plt
import ps_f
import numpy as np
from detect_peaks import detect_peaks
import math

# # Using code from different author, Marcos Duarte, detect_peaks.
print('Type t to run in test-mode (sine-functions), r for regular run.')
inpt0 = input()


if inpt0 == 'r':
    # # Data loading, powerspectrum and corresponding stellar cflux
    print('Which spectrum to use? Insert datafile, excluding .dat')
    inpt1 = input()

    freq = []
    p = []
    alpha = []
    beta = []
    data = ps_f.preader(inpt1 + '.dat')
    for row in data:
        freq.append(row[0])
        p.append(row[1])
        alpha.append(row[2])
        beta.append(row[3])

    print('Which data to use? Use the same data as the spectrum is made from, exclude .dat')
    inpt2 = input()

    time = []
    cflux = []
    data2 = ps_f.preader(inpt2 + '.dat')
    for row in data2:
        time.append(row[0])
        cflux.append(row[1])
    time0 = time[0]
    time[:] = [x-time0 for x in time]  # Correcting time to start at 0


elif inpt0 == 't':  # Using combination of functions: f = 3*sin(2*t) + cos(3*t) - sin(0.5*t)
    steps = 1000
    maxtime = 1000
    time = np.linspace(0, maxtime, steps)

    func = [np.sin(x) for x in time]

    dt = maxtime/steps
    test_nyq = 1 / (2 * dt)
    test_nyq = 2 * np.pi * test_nyq
    resolution = 0.001 * test_nyq / (2*np.pi)
    ps_result = ps_f.frequency_examiner(func, time, [0.159], 0.01, 0.0001)
    ps_result = ps_result[0]
    freq = ps_result[0]
    p = ps_result[1]
    print(len(p))
    alpha = ps_result[2]
    beta = ps_result[3]

    cflux = func
    plt.plot(freq, p)
    plt.show()
else:
    print('Unknown input.')

# # Peak detection
# min peak height perhaps 90?

indices_p = detect_peaks(p, show=True)

# #
signal = [None] * len(indices_p)
signal_tot = [0] * len(time)
for k in range(0, len(indices_p)):
    index = indices_p[k]
    signal[k] = [alpha[index] * np.cos(2*np.pi*freq[k] * x) + beta[index] * np.sin(2*np.pi*freq[k] * x) for x in time]
    signal_temp = signal[k]
    for h in range(0, len(signal_tot)):
        signal_tot[h] = signal_tot[h] + signal_temp[h]

    plt.plot(time, signal[k])

#plt.xlim(0, 100*14000)
plt.show()


plt.plot(time, cflux, 'r')
plt.plot(time, signal_tot, 'b')

#plt.xlim(0, 100*14000)
plt.show()
