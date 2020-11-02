"""
Function "repository" with general functions needed for calculation of power-spectrum.
Separated from run-script ps_run.py for convenience.
"""

import numpy as np
import math
import os
import time as tm
from detect_peaks import detect_peaks
import matplotlib.pyplot as plt
import create_pspectrum


def writer(fname, *args):
    fname = os.path.join('Datafiles', fname + '.dat')
    np.savetxt(fname, np.c_[args])


def reader(fname):
    fname = os.path.join('Datafiles', fname + '.dat')
    arr = np.loadtxt(fname)
    return np.array_split(arr, len(arr[0, :]), axis=1)


def optimal_resolution(t, initial_resolution, sfreq, half_width=None, chunk_size=100):
    """
    Finds a more optimal frequency resolution to reduce the amount of noise introduced
    in the spectrum. Does it by simulating a sine function sampled similarly to the data set
    and examines how the spectral response of it is at different power spectrum resolutions.
    Assumes that the sine frequency used is not very significant for the result. As a sine
    frequency, the frequency of maximum power could be used.
    """
    if half_width is None:
        half_width = sfreq - sfreq/100
    yt = np.sin(2*math.pi*sfreq*t)
    # Initial
    #pres = create_pspectrum.numpy(yt, t, freq_centre=[sfreq], half_width=half_width, resolution=initial_resolution,
    #                              chunk_size=chunk_size)
    pres = create_pspectrum.nufftpy(yt, t, half_width=half_width, resolution=initial_resolution)
    freq_i = pres[0]
    power_i = pres[1]
    # Find area under the curve (sum the power spectrum)
    area_i = np.sum(power_i)
    best_area = area_i
    best_res = initial_resolution
    best_freq = freq_i
    best_p = power_i

    # # # Loop to examine # # # #
    resolution = np.linspace(0.90*initial_resolution, 2.5*initial_resolution, 10)
    for current_res in resolution:
        pres = create_pspectrum.numpy(yt, t, freq_centre=[sfreq], half_width=half_width, resolution=current_res,
                                      chunk_size=chunk_size)[0]
        #  pres = create_pspectrum.nufftpy(yt, t, half_width=half_width, resolution=current_res)
        current_power = pres[1]
        current_area = np.sum(current_power)
        if current_area < best_area:
            best_area = current_area
            best_res = current_res
            best_freq = pres[0]
            best_p = pres[1]

    print('Best resolution: ', best_res)
    print('Initial resolution: ', initial_resolution)
    plt.plot(freq_i/0.000001, power_i, '--')
    plt.plot(best_freq/0.000001, best_p)
    plt.show()
    return best_res



