"""
Function "repository" with general functions needed for calculation of power-spectrum.
Separated from run-script ps_run.py for convenience.
"""

import numpy as np
import math
import os
import time as tm
import nufftpy
from detect_peaks import detect_peaks
import matplotlib.pyplot as plt


def writer(fname, *args):
    fname = os.path.join('Datafiles', fname + '.dat')
    np.savetxt(fname, np.c_[args])


def reader(fname):
    fname = os.path.join('Datafiles', fname + '.dat')
    arr = np.loadtxt(fname)
    return np.array_split(arr, len(arr[0, :]), axis=1)


def create_pspectrum(y, t, freq_centre, half_width, resolution):
    pi = math.pi
    # # Prepare input for use

    # Convert input to numpy array (if not already)
    y = np.asarray(y)
    t = np.asarray(t)

    # Change to angular frequencies (assumes input is in cyclic frequencies)
    freq_centre[:] = [x * (2 * pi) for x in freq_centre]
    half_width = half_width * 2 * pi
    resolution = resolution * 2 * pi

    # Data mean subtraction, to reduce potentially constant elements
    y_mean = np.mean(y)
    y = y - y_mean
    t = t - np.min(t)  # To get first element in t to be 0

    # # Prepare for for-loop

    # Preload results list
    results = [None] * len(freq_centre)

    # Set amount of steps (might be subject to rounding error)
    step_amnt = int((2 * half_width) / resolution)

    # # Run for loop for all frequencies indicated by length of freq_centre

    for k in range(0, len(freq_centre)):
        # Get current frequency centre
        freq_centre_current = freq_centre[k]

        # Create frequency steps
        freq = np.linspace(freq_centre_current - half_width, freq_centre_current + half_width, step_amnt)

        # Reshape freq and t in order to do matrix multiplication
        freq = np.reshape(freq, (len(freq), 1))
        t = np.reshape(t, (1, len(t)))

        # # Zhu Li, do the thing

        # Calculate sine and cosine function part values
        print(freq.shape)
        print(t.shape)
        try:
            product = np.matmul(freq, t)
            sin = np.sin(product)
            cos = np.cos(product)
            sin2 = sin ** 2
            cos2 = cos ** 2
            sincos = sin * cos
            print('sin.shape', sin.shape)

            # Sum together values for different t
            sin = np.sum(y * sin, 1)
            cos = np.sum(y * cos, 1)
            sin2 = np.sum(sin2, 1)
            cos2 = np.sum(cos2, 1)
            sincos = np.sum(sincos, 1)
            print('sin.shape', sin.shape)
            print('Regular run')

        except MemoryError:
            sin = np.zeros(step_amnt)
            cos = np.zeros(step_amnt)
            sin2 = np.zeros(step_amnt)
            cos2 = np.zeros(step_amnt)
            sincos = np.zeros(step_amnt)

            chunk_size = 1000
            freq = np.ascontiguousarray(freq)
            t = np.ascontiguousarray(t)
            for i in range(0, step_amnt, chunk_size):
                end = i + chunk_size
                if end >= step_amnt:
                    end = step_amnt

                product = np.dot(freq[i:end], t)
                # The time-heavy calculation
                sin_temp = np.sin(product)
                cos_temp = np.cos(product)

                sin[i:end] = np.sum(y * sin_temp, 1)
                cos[i:end] = np.sum(y * cos_temp, 1)
                sin2[i:end] = np.sum(sin_temp ** 2, 1)
                cos2[i:end] = np.sum(cos_temp ** 2, 1)
                sincos[i:end] = np.sum(sin_temp * cos_temp, 1)

        # # Calculate alpha and beta components of spectrum, and from them, power of spectrum
        alpha = (sin * cos2 - cos * sincos) / (sin2 * cos2 - np.power(sincos, 2))
        beta = (cos * sin2 - sin * sincos) / (sin2 * cos2 - np.power(sincos, 2))

        power = np.power(alpha, 2) + np.power(beta, 2)

        # # Last loop steps

        # Convert frequency back to cyclic units
        freq = freq / (2 * pi)

        # Save data in results
        results[k] = [freq, power, alpha, beta]

    return results


def clean_procedure(t, y, n_iter, halfwidth, resolution, window=None, mph=1):
    # # Preparation
    # Preload lists for storing found peaks
    p_freq = []
    p_power = []
    p_alpha = []
    p_beta = []
    # Get amount of frequency steps for nufft1
    steps = int((2 * halfwidth) / resolution)
    # Make carbon copy of signal for manipulation
    y_copy = np.copy(y)

    # # Loop
    for i in range(0, n_iter):
        t1 = tm.time()

        # Get cyclic frequencies for nufft1
        freq = nufftpy.nufftfreqs(steps, df=resolution)
        freq = freq[len(freq) // 2:-1]  # Take only positive frequencies
        # Call nufft1 to perform nufft calculation
        harmonic_content = nufftpy.nufft1(t, y_copy, steps, df=(resolution * 2 * np.pi))
        harmonic_content = harmonic_content[len(harmonic_content)//2:-1]
        # Create window (if None) and apply to spectrum to only observe relevant area
        if window is None:
            window = range(0, len(harmonic_content))
        harmonic_content = harmonic_content[window]
        freq = freq[window]
        # Calculate power
        spectral_power = harmonic_content.real ** 2 + harmonic_content.imag ** 2

        # Detect peaks in power spectrum, and find power, frequency, alpha and beta values
        peaks = detect_peaks(spectral_power, mph=mph)
        peaks_power = spectral_power[peaks]
        peaks_alpha = harmonic_content.real[peaks]
        peaks_beta = harmonic_content.imag[peaks]
        peaks_freq = freq[peaks]

        # Find highest peak
        max_indx = np.argmax(peaks_power)
        max_freq = peaks_freq[max_indx]
        max_power = peaks_power[max_indx]
        max_alpha = peaks_alpha[max_indx]
        max_beta = peaks_beta[max_indx]

        # plt.plot(freq, spectral_power)
        # plt.plot(peaks_freq, peaks_power, 'r*', markersize=4)
        # plt.plot(max_freq, max_power, 'b*', markersize=6
        # plt.show()

        # Calculate harmonic signal corresponding to highest peak in power spectrum
        max_signal = (max_alpha * np.cos(2*np.pi*max_freq * t) + max_beta * np.sin(2*np.pi*max_freq * t)) * 2

        # plt.plot(t, y_copy, linewidth=0.5)
        # plt.plot(t, max_signal, '--', linewidth=0.3)
        # plt.show()

        # Subtract calculated signal from y_copy and save the peak used
        y_copy -= max_signal
        p_power.append(max_power)
        p_alpha.append(max_alpha)
        p_beta.append(max_beta)
        p_freq.append(max_freq)

        t2 = tm.time()
        print(t2-t1)

    # # Result
    # Get cyclic frequencies for nufft1
    freq = nufftpy.nufftfreqs(steps, df=resolution)
    freq = freq[len(freq) // 2:-1]
    # Call nufft1 to perform nufft calculation for unchanged signal y
    harmonic_content = nufftpy.nufft1(t, y, steps, df=(resolution * 2 * np.pi))
    harmonic_content = harmonic_content[len(harmonic_content) // 2:-1]
    # Calculate spectral power
    spectral_power = harmonic_content.real ** 2 + harmonic_content.imag ** 2

    # Plot power spectrum and peaks found by cleaning
    plt.plot(freq, spectral_power)
    plt.plot(p_freq, p_power, 'r*', markersize=4)
    plt.show()

    # Save peaks found by cleaning in .dat file
    writer('clean_peaks', p_freq, p_power)

    return p_freq, p_power, p_alpha, p_beta
