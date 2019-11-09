"""
Script to generate power-spectra by using light-curve data from data_filterer.py. Calls functions from the function
'repository' ps_f.py to do this.

Summary:
-   Loads data from .dat file.
-   Either:
           - Makes a crude power-spectrum calculation by calling psspectrum (3000 frequencies), and saves in data file
           - Loads previously created crude power spectrum from .dat file
-   Plots crude power spectrum, and asks for graphical input to select examination frequencies for function
    ps_f.frequency_examiner
-   Generates one or more high-resolution power spectra by using ps_f.frequency, with a resolution of 0.01 microHz
    (and a width of 1000 microHz).
-   Saves generated power-spectra to data files, with the original filename + '_ps_ex' + str(i), where i is the number
    of power spectrum saved this run (starting from 0).
"""
# # Run file to create power spectrum. Calls functions from ps_f.py. Currently analyzing data from nu indi/indus.

import math
import ps_f
import matplotlib.pyplot as plt

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


# # Input to decide whether to skip calculation and load previously made power spectrum
print('Do you want to reload saved data? (y/n) If not, new power spectrum will be made.')
inpt2 = input()
if inpt2 == "y":
    repeat = False
elif inpt2 == "n":
    repeat = True
else:
    print('Did not enter correct letter. Assuming you meant y.')
    repeat = False


if repeat is True:  # Makes a 'short' calculation of the power spectrum, in order to show a preview
    # TESS data info
    dt = 20  # seconds, TESS sample rate
    nyquist = 1/(2*dt)  # Nyquist frequency
    angular_nyq = 2*math.pi*nyquist  # Angular nyquist frequency

    steps = 3000  # Amount of frequencies to try

    # # Function call, crude power spectrum made
    k = 0.9  # max frequency factor (k*nuqvist freq)
    results = ps_f.pspectrum(cflux, tid, angular_nyq, steps, k)  # gives frequency back in cyclic units
    results = results[0]
    [freq, p, alpha, beta] = [results[0], results[1], results[2], results[3]]  # p=power spectrum, alpha beta is parts

elif repeat is False:
    data = ps_f.preader('ps_crude.dat')  # Read function call, reads the file 'ps_crude.dat' for an existing spectrum

    # # Preload variables
    freq = []
    p = []
    alpha = []
    beta = []

    # # Append data to variables
    for row in data:
        freq.append(row[0])
        p.append(row[1])
        alpha.append(row[2])
        beta.append(row[3])


# # Plotting of crude power spectrum
plt.plot(freq, p)
plt.xlabel('Cyclical frequency [Hz]')
plt.ylabel('Power Spectrum')
plt.show()


# # Saving to file, only if new data is made
if repeat is True:
    ps_f.pwriter('ps_crude', freq, p, alpha, beta)


# # Frequency Conversion
muHz = 0.000001  # Variable for converting to microHz
freq_mHz = [x / muHz for x in freq]  # Frequencies converted to microHz


# # Power spectrum plot and graphical selection of point to examine around
plt.plot(freq_mHz, p)
plt.xlabel('Cyclical frequency [muHz]')
plt.ylabel('Power Spectrum')
plt.title('Select frequencies to examine')
coords = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)  # coordinates of cursor
plt.close()

# # Convert input(s) to list elements indicating only the frequency selected
frequencies = [None] * len(coords)  # Preload list of frequencies to examine around
for k in range(0, len(coords)):
    coord_select = coords[k]
    frequencies[k] = coord_select[0] * muHz  # Frequency to examine around (necessary for function call)


# # Frequency examination
resolution = 0.01 * muHz  # 0.01 normal
# dt = 120  # seconds, TESS sample rate
# nyquist = 1/(2*dt)  # Nyquist frequency
# resolution = nyquist / 10000
examined_results = ps_f.create_pspectrum(cflux, tid, frequencies, 2000 * muHz, resolution)  # Function call to
# examine around the indicated points. Is done in a broad interval (2*500 mikroHz), so the practically all of the
# relevant parts of the interval. Generally only 1 frequency to examine around is used because of this.

plt.plot(freq_mHz, p, 'r')  # plot the crude power spectrum to compare with below plots

for h in range(0, len(examined_results)):  # loop to plot and save all the spectrum intervals examined
    # # Pull current loop results from list
    temp_results = examined_results[h]
    [freq_exam, p_exam, alp_exam, bet_exam] = [temp_results[0], temp_results[1], temp_results[2], temp_results[3]]
    freq_exam_mHz = [x / muHz for x in freq_exam]

    # # Plot current loop results, and save them to file with added number indicating current loop
    plt.plot(freq_exam_mHz, p_exam)
    ps_f.pwriter(inpt1+'_ps_ex'+str(h), freq_exam, p_exam, alp_exam, bet_exam)

plt.xlabel('Cyclical frequency (muHz)')
plt.ylabel('Power Spectrum')
plt.show()
