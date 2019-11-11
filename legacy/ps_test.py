"""
Script used to test the functions used to generate power spectra, ps_f.pspectrum and ps_f.frequency_examiner.
It does this by feeding them a regular sine function. In turn, the theoretical fourier transform is known for this
sine function and, as such, an expected power spectrum is known. This is then compared graphically with the results from
the two functions, in order to verify their validity.
"""

import ps_f
import numpy as np
import math
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 22})


xtest = np.linspace(0, 1000, 1000)
ytest = [np.sin(x) for x in xtest]

test_dx = 1
test_nyq = 1/(2*test_dx)
test_nyq = 2*math.pi*test_nyq

spacing = 3000

result = ps_f.pspectrum(ytest, xtest, test_nyq, spacing, 0.5)
result = result[0]

freq = result[0]
p = result[1]



# # Plotting
ang_freq = [x*2*math.pi for x in freq]
# # Comparison theoretical power spectrum
p_theory = [4*np.sin(xtest[-1]/2 * (x - 1))**2 / (xtest[-1]**2 * (x-1)**2) for x in ang_freq]  # normalized

plt.plot(ang_freq, p, 'r')
plt.plot(ang_freq, p_theory, 'b--')
plt.xlabel('Angular frequency [rad*Hz]')
plt.ylabel('Power Spectrum')
plt.legend(['Calculated power-spectrum', 'Normalized theoretical power-spectrum'])
plt.show()



# # Frequency selection
plt.plot(freq, p)
plt.xlabel('Cyclical frequency [Hz]')
plt.ylabel('Power Spectrum')
coords = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)


frequencies = [None] * len(coords)
for k in range(0, len(coords)):
    coord_select = coords[k]
    frequencies[k] = coord_select[0]


# # Frequency examination
examined_results = ps_f.frequency_examiner(ytest, xtest, frequencies, 0.01, 0.00001)
plt.plot(freq, p, 'r')
for h in range(0, len(examined_results)):
    temp_results = examined_results[h]
    plt.plot(temp_results[0], temp_results[1], 'b')

plt.show()

plt.plot(ang_freq, p, 'r')
first_res = examined_results[0]
first_freq = first_res[0]
first_p = first_res[1]
first_freq_ang = [x * 2 * math.pi for x in first_freq]

plt.plot(first_freq_ang, first_p, 'b--')
plt.legend(['Rough power-spectrum', 'Detailed power-spectrum'])
plt.xlim([0.95, 1.05])
plt.xlabel('Angular frequency [rad*Hz]')
plt.ylabel('Power Spectrum')
plt.show()


plt.plot(first_freq_ang, first_p, 'r')
first_res = examined_results[0]
first_freq = first_res[0]
first_p = first_res[1]

first_freq_ang = [x * 2 * math.pi for x in first_freq]
first_p_theory = [4*np.sin(xtest[-1]/2 * (x - 1))**2 / (xtest[-1]**2 * (x-1)**2) for x in first_freq_ang]  # normalized

plt.plot(first_freq_ang, first_p_theory, 'b--')
plt.legend(['Calculated power-spectrum', 'Normalized theoretical power-spectrum'], loc=2)
plt.xlim([0.94, 1.05])
plt.xlabel('Angular frequency [rad*Hz]')
plt.ylabel('Power Spectrum')
plt.show()
