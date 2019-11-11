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


def reader_legacy(filename):
    """
    Function to read .dat files. Assumes that there are 4 columns in the file. Generally used to read lightcurve data
    sets generated by "data_filterer.py" (time, cflux, median cflux, variance)

    Returns:
        a:  list of lists, where each list has 4 values, one matching the different variables contained in the .dat-file
            can be seen as 4 columns, where the number of lists are number of rows in a matrix.

    Variables:
        filename:   string with the name of the file, including .dat extension
    """

    a = []
    f = open(os.path.join('Datafiles', filename), 'r')
    for line in f:
        columns = line.split()
        columns[0] = float(columns[0])
        columns[1] = float(columns[1])
        columns[2] = float(columns[2])
        columns[3] = float(columns[3])
        a.append(columns)
    f.close()
    return a


def writer(fname, *args):
    fname = os.path.join('Datafiles', fname + '.dat')
    np.savetxt(fname, np.c_[args])


def reader(fname):
    fname = os.path.join('Datafiles', fname + '.dat')
    arr = np.loadtxt(fname)
    return np.array_split(arr, len(arr[0, :]), axis=1)


def pfunctions(func, t, max_freq, spacing, k):
    """
    Function to calculate least-sine squared sums in order to generate power-spectrum (calculates dft sums for
    different, equally spaced frequencies with possible data gaps).
    Generally used to calculate a crude and quick power-spectrum, in order to decide where it is relevant to improve
    resolution with the function frequency_examiner.

    Returns:
        res: list of lists, res = [nu_res, s_res, c_res, ss_res, cc_res, sc_res], contains a list with each of the sum
             values for different frequencies
             nu_res = list of frequencies
             s_res = list of s sums, one for each frequency
             c_res = similar
             ...

             For a description of the sums, see
                S. Frandsen et al., Calculation of the power spectrum, Astronomy & Astrophysics, 311:123-134, 1995,
                Section 4.2 of article, pages 125-126

    Variables:
        func: List of floats, lightcurve data, y-values, normally corrected flux PPM.
        t: List of floats, time data for func, x-values, in seconds.
        max_freq: Float, the highest frequency the function should generate sums for. Often the nyquist frequency.
        spacing: Integer, the desired amount of steps (related to frequency resolution).
        k: Float, a factor to multiply the max_freq with. To be used f.ex. when max_freq is set to the nyquist
           frequency, but you only want to generate sums to 0.5*nyquist (then k=0.5).
    """
    t[:] = [x - min(t) for x in t]  # to get first element to be 0
    d_nu = k * max_freq / spacing
    nu = 5 * d_nu  # Skips over first 4 steps, as the frequency will be too low (near nyquist)

    # # Preload variables as lists, to be used for sum functions (one list index for each frequency)
    s_res = [None] * spacing
    c_res = [None] * spacing
    ss_res = [None] * spacing
    cc_res = [None] * spacing
    sc_res = [None] * spacing
    nu_res = [None] * spacing

    # # Set starting list index w
    w = 0

    while nu <= (max_freq * k):  # Loop over frequencies, stops when nu has reached max_freq (often nyquist freq) times
        # a factor k

        # # Preload lists for function values
        s_list = [0] * len(t)
        c_list = [0] * len(t)
        ss_list = [0] * len(t)
        cc_list = [0] * len(t)
        sc_list = [0] * len(t)

        for i in range(0, len(t)):  # Loops over elements in t and creates values for functions (previous lists)
            sin = np.sin(nu * t[i])
            cos = np.cos(nu * t[i])
            s_list[i] = func[i] * sin
            c_list[i] = func[i] * cos
            ss_list[i] = sin**2
            cc_list[i] = cos**2
            sc_list[i] = sin * cos

        # # Sums together all function values for a single frequency, adds to result lists
        s_res[w] = np.sum(s_list)
        c_res[w] = np.sum(c_list)
        ss_res[w] = np.sum(ss_list)
        cc_res[w] = np.sum(cc_list)
        sc_res[w] = np.sum(sc_list)
        nu_res[w] = nu  # Saves current frequency

        # # Advance to next frequency and index
        nu += d_nu
        w += 1

    m = len(s_res)  # Length of lists, to be used to correct below

    for k in range(1, m+1):  # Float/integer problem, not using integer step size d_nu, deleted additional empty
        # list elements in loop going from end to start
        if ss_res[m-k] is None:
            del s_res[m-k]
            del c_res[m-k]
            del ss_res[m-k]
            del cc_res[m-k]
            del sc_res[m-k]
            del nu_res[m-k]

    res = [nu_res, s_res, c_res, ss_res, cc_res, sc_res]  # All results from function
    return res


def psummer(s, c, ss, cc, sc):
    """
    A function to multiply the least-sine squared sums together in the correct way, in order to generate the power
    spectrum, and the alpha and beta values. To be used in conjunction with the function pfunctions (which calculates
    the sum functions) or similar.

    Returns:
        res: list of lists, [p, alpha, beta]
             p: list of floats, containing power spectrum values, each with a matching frequency known from the sum
                functions
             alpha: list of floats, the alpha part of the signal
             beta: list of floats, the beta signal part. alpha**2 + beta**2 = p

    Variables:
        s, c, ss, cc, sc: lists of floats, containing sum function values. See documentation for the function pfunctions
                          descriptions of them.
    """
    # # Preload lists
    alpha = [0] * len(s)
    beta = [0] * len(s)
    p = [0] * len(s)

    for x in range(0, len(alpha)):  # Loop to create power-spectrum and components
        alpha[x] = (s[x] * cc[x] - c[x] * sc[x]) / (ss[x] * cc[x] - sc[x] ** 2)
        beta[x] = (c[x] * ss[x] - s[x] * sc[x]) / (ss[x] * cc[x] - sc[x] ** 2)
        p[x] = alpha[x]**2 + beta[x]**2

    res = [p, alpha, beta]
    return res


def pspectrum(func, t, max_freq, steps, max_freq_factor):
    """
    Function that prepares data for use with pfunctions and psummer, and calls both to produce a power spectrum, then
    converts the result to cyclic frequencies.

    Returns:
        res: list of lists, [res1, res2]
             res1: list of lists, [freq, p, alpha, beta]
                   freq: cyclical frequency for power spectrum
                   p, alpha, beta: power spectrum values, and its components, see psummer for description.
             res2: the returns from pfunctions, see pfunctions for description

    Variables:
        func, t, max_freq, steps, max_freq_factor (k): same as for pfunctions, see pfunctions for description.

    """
    # Function to call the two functions pfunctions, psummer, and generate a power spectrum from it in cyclic freq.

    # # Data mean subtraction, to reduce additional constant element (frequency 0)
    func_mean = np.mean(func)
    func[:] = [x - func_mean for x in func]

    # # pfunctions call to create sum functions
    pf = pfunctions(func, t, max_freq, steps, max_freq_factor)

    # # psummer call to create power-spectrum using sum functions
    ps = psummer(pf[1], pf[2], pf[3], pf[4], pf[5])

    # # Retrieve angular frequency, ps and components as lists
    [nu, p, alpha, beta] = [pf[0], ps[0], ps[1], ps[2]]

    # # Convert to cyclic frequency
    freq = [x / (2*math.pi) for x in nu]

    # # Group results
    res1 = [freq, p, alpha, beta]
    res2 = pf
    res = [res1, res2]

    return res


def frequency_examiner(func, t, frequencies, halfwidth, resolution):
    """
    Function with the same overall purpose as the pspectrum function. It calculates power spectrum sum functions (see
    pfunctions for reference), and uses this to calculate the power spectrum (calls psummer). However, this function
    is specified to generate power spectra around one or more indicated frequencies, with an indicated half-width, and a
    resolution instead of number of steps. In this sense, it is more intuitive, and more flexible, as it can be used
    to examine only a specific, smaller region, but also the whole region like the others. It is also easier to
    indicate a specific cyclic frequency resolution to be used, instead of indicating an amount of steps (where the
    resolution then depends on the interval size). This, however, also means that there is less control on the amount
    of time that the function should spend.

    Returns:
        results: list of lists of lists, 1 list for each element in the frequencies variable.
                 Each list: [freq_res, p, alpha, beta], data values for the region around the examined frequency
                            freq_res: List of floats, power spectrum frequencies in cyclical units
                            p: list of floats, power spectrum values for the frequencies of freq_res
                            alpha, beta: lists of floats, signal components for the frequencies of freq_res

    Variables:
        func: List of floats, lightcurve data, y-values, normally corrected flux PPM.
        t: List of floats, time data for func, x-values, in seconds.
        frequencies: List of floats, detailing specific frequencies to examine around. Usually, only 1 is used, and it
                     is usually chosen by first making a crude power spectrum with pspectrum function. Cyclical units.
        halfwidth:   Float, cyclical frequency difference between the lowest frequency, and the frequency specified by
                     'frequencies'. It is half the width of the frequency region examined.
        resolution:  Float, frequency resolution in cyclical frequency units (f.ex. 0.01Hz).
    """

    # # Change to angular frequencies (assumes input is in cyclic frequencies)
    frequencies[:] = [x * (2*math.pi) for x in frequencies]
    halfwidth = halfwidth * 2 * math.pi
    resolution = resolution * 2 * math.pi

    # # Data mean subtraction, to reduce additional constant element (frequency 0)
    func_mean = np.mean(func)
    func[:] = [x - func_mean for x in func]

    t[:] = [x - min(t) for x in t]  # to get first element to be 0

    # # Total amount of steps (might be subject to a rounding error)
    steps = int((2 * halfwidth) / resolution)

    # # Preload list for results
    results = [None] * len(frequencies)

    # # Full examination loop to do power-spectrum calculation for all indicated frequencies (can be a list)
    for k in range(0, len(frequencies)):

        # # Pull specific frequency to examine around
        freq = frequencies[k]

        # # Set starting frequency for while loop
        nu = freq - halfwidth

        # # Set starting list index for while loop
        w = 0

        # # Preload resulting sum function lists, with multiple "extra" elements to assure index length is not exceeded
        s_res = [None] * (steps+5)  # We sort out extra list elements in end of loop
        c_res = [None] * (steps+5)
        ss_res = [None] * (steps+5)
        cc_res = [None] * (steps+5)
        sc_res = [None] * (steps+5)
        nu_res = [None] * (steps+5)

        # # Loop over frequencies, stops when nu has passed freq by + 1halfwidth
        while nu <= freq + halfwidth:

            # # Preload lists for function values for current loop frequency
            s_list = [0] * len(t)
            c_list = [0] * len(t)
            ss_list = [0] * len(t)
            cc_list = [0] * len(t)
            sc_list = [0] * len(t)

            # # Loop to calculate function values for the whole time-set (for current loop frequency)
            for i in range(0, len(t)):
                sin = np.sin(nu * t[i])
                cos = np.cos(nu * t[i])
                s_list[i] = func[i] * sin
                c_list[i] = func[i] * cos
                ss_list[i] = sin ** 2
                cc_list[i] = cos ** 2
                sc_list[i] = sin * cos

            # # Sum together all function values for current frequency, adds to result lists
            s_res[w] = np.sum(s_list)
            c_res[w] = np.sum(c_list)
            ss_res[w] = np.sum(ss_list)
            cc_res[w] = np.sum(cc_list)
            sc_res[w] = np.sum(sc_list)
            nu_res[w] = nu  # Saves current frequency

            # # Advance to next frequency and index
            nu += resolution
            w += 1

        # Float/integer problem, not using integer step size d_nu, deleting additional empty list elements in loop
        # going from end to start
        m = len(s_res)
        for j in range(1, m + 1):
            if ss_res[m - j] is None:
                del s_res[m - j]
                del c_res[m - j]
                del ss_res[m - j]
                del cc_res[m - j]
                del sc_res[m - j]
                del nu_res[m - j]

        # # Group all results from while loop (sum functions)
        pf = [nu_res, s_res, c_res, ss_res, cc_res, sc_res]

        # # Call psummer to create power-spectrum from sum functions
        ps = psummer(pf[1], pf[2], pf[3], pf[4], pf[5])

        # # Pull frequency, power-spectrum, and components from psummer list
        [freq_res, p, alpha, beta] = [pf[0], ps[0], ps[1], ps[2]]

        # # Convert frequency back to cyclic units
        freq_res[:] = [x / (2*math.pi) for x in freq_res]

        # # Save to results list
        results[k] = [freq_res, p, alpha, beta]

    return results


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

        # plt.plot(freq, spectral_power)
        # plt.plot(peaks_freq, peaks_power, 'r*', markersize=4)
        # plt.show()

        # Find highest peak
        max_indx = np.argmax(peaks_power)
        max_freq = peaks_freq[max_indx]
        max_power = peaks_power[max_indx]
        max_alpha = peaks_alpha[max_indx]
        max_beta = peaks_beta[max_indx]

        # Calculate harmonic signal corresponding to highest peak in power spectrum
        max_signal = (max_alpha * np.cos(2*np.pi*max_freq * t) + max_beta * np.sin(2*np.pi*max_freq * t)) * 2

        # plt.plot(t, y_copy, linewidth=0.5)
        # plt.plot(t, max_signal, '--', linewidth=0.3)
        # plt.show()

        # Subtract calculated signal from y_copy and save the peak used
        y_copy -= max_signal
        p_power.append(max_power)
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


def pwriter(filename, freq, p, alpha, beta):
    # # Used to write power spectrum data to a file. Assumes freq, p, alpha, beta are lists of equal length.

    # # Clear file to be used (if it exists, otherwise it creates an empty file by the name)
    open(os.path.join('Datafiles', filename + '.dat'), 'w').close()

    # # Loop over all elements, save matching elements to a single line (freq[i], p[i], alpha[i], beta[i])
    k = len(freq)
    for i in range(0, k):
        content = str(freq[i]) + '   ' + str(p[i]) + '   ' + str(alpha[i]) + '   ' + str(beta[i]) + '\n'
        with open(os.path.join('Datafiles', filename + '.dat'), 'a') as f:
            f.write(content)


def writer2col(filename, col1, col2):
    # # Used to write 2 lists of same size to a file. Expects list of floats, strings or integers.

    # # Clear file to be used (if it exists, otherwise it creates an empty file by the name)
    open(os.path.join('Datafiles', filename + '.dat'), 'w').close()

    # # Loop over all elements, save matching elements to a single line (freq[i], p[i], alpha[i], beta[i])
    k = len(col1)
    for i in range(0, k):
        content = str(col1[i]) + '   ' + str(col2[i]) + '\n'
        with open(os.path.join('Datafiles', filename + '.dat'), 'a') as f:
            f.write(content)


def preader(filename):
    A = []
    f = open(os.path.join('Datafiles', filename), 'r')
    for line in f:
        columns = line.split()
        columns[0] = float(columns[0])
        columns[1] = float(columns[1])
        columns[2] = float(columns[2])
        columns[3] = float(columns[3])
        A.append(columns)
    f.close()
    return A


def dreader(filename):
    A = []
    f = open(os.path.join('Datafiles', filename), 'r')
    for line in f:
        columns = line.split()
        columns[0] = float(columns[0])
        columns[1] = float(columns[1])
        A.append(columns)
    f.close()
    return A


def gauss(x, amp, cen, wid, const):
    return amp * np.exp(-(x-cen)**2 / wid) + const


def ngauss(x, a, b, c, d, n):  # a, b, c must be lists with n elements
    result = 0
    for i in range(0, n):
        result += gauss(x, a[i], b[i], c[i], d[i])
    return result


def neg_exp(x, a, b, c):
    return a * np.exp(-b * x) + c


def lorentz(x, width, center, white):
    return (1/np.pi) * (0.5 * width) / ((x-center)**2 + (0.5*width)**2) + white


def lorentz2(x, amplitude, sigma, center, c):
    return c + (amplitude / np.pi) * (sigma/((x - center)**2 + sigma**2))


def lorentz3(x, amplitude, sigma, center):
    return (amplitude / np.pi) * (sigma/((x - center)**2 + sigma**2))


def nlorentz(x, width, center, white, n):
    result = 0
    for i in range(0, n):
        result += lorentz(x, width[i], center[i], white[i])
    return result
