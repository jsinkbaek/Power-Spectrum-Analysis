import ps_f
import time as tm
from detect_peaks import detect_peaks
import numpy as np
import matplotlib.pyplot as plt
import create_pspectrum


def nufft(t, y, n_iter, halfwidth, resolution, window=None, mph=1):
    """
    Perform CLEAN procedure using FFT (NUFFT)
    """
    # # Preparation
    import nufftpy
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
        try:
            max_indx = np.argmax(peaks_power)
            max_freq = peaks_freq[max_indx]
            max_power = peaks_power[max_indx]
            max_alpha = peaks_alpha[max_indx]
            max_beta = peaks_beta[max_indx]
        except ValueError:
            break
        # plt.figure()
        # plt.plot(freq, spectral_power)
        # plt.plot(peaks_freq, peaks_power, 'r*', markersize=4)
        # plt.plot(max_freq, max_power, 'b*', markersize=6)
        # plt.show(block=False)

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
    plt.figure()
    plt.plot(freq, spectral_power)
    plt.plot(p_freq, p_power, 'r*', markersize=4)
    plt.show()

    # Save peaks found by cleaning in .dat file
    ps_f.writer('clean_peaks', p_freq, p_power)

    return p_freq, p_power, p_alpha, p_beta


def numpy(t, y, n_iter, freq_centre, fft_half_width, resolution, np_half_width=10, chunk_size=100, window=None, mph=1):
    """
    CLEAN procedure using numpy matrix multiplication on a small frequency range, using nufftpy to find the approximate
    frequency of the highest peak beforehand.
    """
    # # # Preparation # # #
    import nufftpy
    p_freq = []
    p_power = []
    p_alpha = []
    p_beta = []
    freq_centre = [freq_centre]
    # Get amount of frequency steps for nufft1
    steps = int((2 * fft_half_width) / resolution)
    # Make carbon copy of signal for manipulation
    y_copy = np.copy(y)

    # # # # CLEAN LOOP # # # #
    for i in range(0, n_iter):
        t1 = tm.time()

        # # # Do FFT to find expected highest peak area # # #
        # Get cyclic frequencies for nufft1
        freq = nufftpy.nufftfreqs(steps, df=resolution)
        freq = freq[len(freq) // 2:-1]  # Take only positive frequencies
        # Call nufft1 to perform nufft calculation
        harmonic_content = nufftpy.nufft1(t, y_copy, steps, df=(resolution * 2 * np.pi))
        harmonic_content = harmonic_content[len(harmonic_content) // 2:-1]
        spectral_power = harmonic_content.real ** 2 + harmonic_content.imag ** 2
        # Detect peaks in power spectrum, and find power, frequency, alpha and beta values
        peaks = detect_peaks(spectral_power, mph=mph)
        peaks_power = spectral_power[peaks]
        peaks_freq = freq[peaks]
        # Find highest peak
        try:
            max_indx = np.argmax(peaks_power)
            max_freq = peaks_freq[max_indx]
        except ValueError:
            break
        t2 = tm.time()
        print('fft time: ', t2-t1)

        # # # Do sine-cosine fitting to estimate power spectrum using numpy matrix multiplication # # #
        spectral_res = create_pspectrum.numpy(y_copy, t, half_width=np_half_width, freq_centre=max_freq,
                                              resolution=resolution, chunk_size=chunk_size)
        spectral_power = spectral_res[1]
        # cut data to window
        if window is None:
            window = range(0, len(spectral_power))
        spectral_power = spectral_power[window]
        # # Do peak finding # #
        peaks = detect_peaks(spectral_power, mph=mph)
        peaks_power = spectral_power[peaks]
        peaks_alpha = spectral_res[2][peaks]
        peaks_beta = spectral_res[3][peaks]
        peaks_freq = spectral_res[0][peaks]
        # Find highest peak
        try:
            max_indx = np.argmax(peaks_power)
            max_freq = peaks_freq[max_indx]
            max_power = peaks_power[max_indx]
            max_alpha = peaks_alpha[max_indx]
            max_beta = peaks_beta[max_indx]
        except ValueError:
            break

        # # # Calculate harmonic signal corresponding to highest peak in power spectrum # # #
        max_signal = (max_alpha * np.cos(2*np.pi*max_freq * t) + max_beta * np.sin(2*np.pi*max_freq * t)) * 2

        # # Subtract calculated signal from y_copy and save the peak used # #
        y_copy -= max_signal
        p_power.append(max_power)
        p_alpha.append(max_alpha)
        p_beta.append(max_beta)
        p_freq.append(max_freq)
        t3 = tm.time()
        print('numpy time: ', t3-t2)

    # # # # Results # # # #
    # Calculate original power spectrum
    spectral_res = create_pspectrum.numpy(y, t, freq_centre=freq_centre, half_width=fft_half_width,
                                          resolution=resolution, chunk_size=chunk_size)
    # Plot spectrum and found peaks
    plt.figure()
    plt.plot(spectral_res[0], spectral_res[1])
    plt.plot(p_freq, p_power, 'r*', markersize=4)
    plt.show()
    # Save peaks found by cleaning in .dat file
    ps_f.writer('clean_peaks_numpy', p_freq, p_power)
    return p_freq, p_power, p_alpha, p_beta


def cuda(t, y, n_iter, freq_centre, half_width, resolution, chunk_size=100, window=None, mph=1):
    """
    CLEAN procedure utilizing create_pspectrum.cuda
    """
    # Preload lists for storing found peaks
    p_freq = []
    p_power = []
    p_alpha = []
    p_beta = []
    freq_centre = [freq_centre]
    # Get amount of frequency steps wanted
    steps = int((2 * half_width) / resolution)
    # Make carbon copy of signal for manipulation
    y_copy = np.copy(y)

    # # # # CLEAN LOOP # # # #
    for i in range(0, n_iter):
        t1 = tm.time()

        # # Perform sine-cosine fit with gpu to get power spectrum # #
        spectral_res = create_pspectrum.cuda(y_copy, t, freq_centre=freq_centre, half_width=half_width,
                                             resolution=resolution, chunk_size=chunk_size)
        spectral_power = spectral_res[1]
        # cut data to window
        if window is None:
            window = range(0, len(spectral_power))
        spectral_power = spectral_power[window]

        # # Do peak finding # #
        peaks = detect_peaks(spectral_power, mph=mph)
        peaks_power = spectral_power[peaks]
        peaks_alpha = spectral_res[2][peaks]
        peaks_beta = spectral_res[3][peaks]
        peaks_freq = spectral_res[0][peaks]

        # Find highest peak
        try:
            max_indx = np.argmax(peaks_power)
            max_freq = peaks_freq[max_indx]
            max_power = peaks_power[max_indx]
            max_alpha = peaks_alpha[max_indx]
            max_beta = peaks_beta[max_indx]
        except ValueError:
            break

        # Calculate harmonic signal corresponding to highest peak in power spectrum
        max_signal = (max_alpha * np.cos(2*np.pi*max_freq * t) + max_beta * np.sin(2*np.pi*max_freq * t)) * 2

        # # Subtract calculated signal from y_copy and save the peak used # #
        y_copy -= max_signal
        p_power.append(max_power)
        p_alpha.append(max_alpha)
        p_beta.append(max_beta)
        p_freq.append(max_freq)

        t2 = tm.time()
        print(t2-t1)

    # # # # Results # # # #
    # Calculate original power spectrum
    spectral_res = create_pspectrum.cuda(y, t, freq_centre=freq_centre, half_width=half_width, resolution=resolution,
                                         chunk_size=chunk_size)
    # Plot spectrum and found peaks
    plt.figure()
    plt.plot(spectral_res[0], spectral_res[1])
    plt.plot(p_freq, p_power, 'r*', markersize=4)
    plt.show()
    # Save peaks found by cleaning in .dat file
    ps_f.writer('clean_peaks_cuda', p_freq, p_power)
    return p_freq, p_power, p_alpha, p_beta
