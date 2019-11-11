# # # # # # #
# This is a function "repository", that stores some functions used for the asterseismology-scripts.
# # # # # # #

# # Import Statements for functions
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import signal
from pandas import Series


def median_noise(x_values, y_values, **kwargs):
    """
    This function calculates median noise levels for one or multiple separate data intervals.

    Variables: x_values, y_values. Must each be lists of floats or integers, lists must be of same length.

    Returns:
        median (list with a median value for each interval)
        intervals (start and end list index of each interval)
        y_intervals (list of lists, each with all y_values in the interval)
        x_intervals (same as y_intervals, just with x_values)

    Optional variables:
        number_of_intervals (can be specified when calling the function. Otherwise, a prompt is issued)

    """

    number_of_intervals = kwargs.get('number_of_intervals', None)  # Used to make it possible to skip first manual input
    if number_of_intervals is None:
        print('How many different median')
        inpt1 = input()
        number_of_intervals = int(inpt1)
    else:
        number_of_intervals = int(number_of_intervals)

    # # Graphical input, indicate which intervals the power spectrum should be split in to calculate separate medians
    plt.plot(y_values)
    plt.title('Select interval start and end for each interval (2*number of interval)')
    indices = plt.ginput(n=2*number_of_intervals, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    plt.close()

    # # Preload list of lists for each interval
    intervals = [[] for i in range(0, number_of_intervals)]
    # # Starting index from graphical input result
    k = 0

    # Loop over number of intervals, to fill the list "intervals"
    for i in range(0, number_of_intervals):
        # # Get start and endpoint for current interval
        [coord1, coord2] = [indices[i+k], indices[i+k+1]]

        # # Only keep x-axis input
        [coord1, coord2] = [coord1[0], coord2[0]]
        if i == 0:
            coord1 = 0
        # # Save interval start and endpoint as integers in interval list
        intervals[i] = [int(coord1), int(coord2)]

        # # Advance k (Ensures the loop skips a beat, so 2*i)
        k += 1

    # # Preload list of lists, each with x, y, or median data for a specific interval
    y_intervals = [[] for i in range(0, number_of_intervals)]
    x_intervals = [[] for i in range(0, number_of_intervals)]
    median = [[] for i in range(0, number_of_intervals)]

    # # Loop to fill lists in the pre-loaded lists with data for each interval
    for i in range(0, number_of_intervals):
        # Select current interval
        interval_select = intervals[i]

        # Fill current lists
        y_intervals[i] = y_values[interval_select[0]:interval_select[1]]
        x_intervals[i] = x_values[interval_select[0]:interval_select[1]]
        median[i] = np.median(y_intervals[i])

    # # Plotting of all the medians, along with their intervals
    fig, ax = plt.subplots()
    ax.set_prop_cycle(plt.cycler('color', ['red', 'blue', 'blue', 'red', 'green', 'magenta']))
    legend_list = []
    for i in range(0, number_of_intervals):
        plt.plot(x_intervals[i], y_intervals[i], '*')
        plt.plot(x_intervals[i], [median[i]] * len(x_intervals[i]))
        legend_list = legend_list + ['Region' + str(i+1), 'Median noise of region' + str(i+1)]
    plt.legend(legend_list)
    plt.xlabel('Frequency [' + r'mu' + 'Hz]')
    plt.ylabel('Power spectrum')
    plt.show()

    return [median, intervals, y_intervals, x_intervals]


def func_over_noise(median, y_intervals):
    """
    Function used to produce signal/noise ratio, given multiple data intervals and a median noise value for each.
    Is designed to use output from the function median_noise as standard.

    Returns:
        list with signal/noise values

    Variables:
        median (list of median values, one for each interval)
        y_intervals (list of lists, containing multiple intervals of datasets. Visual: [[interval1], [interval2], ...].
                     each data-set is expected to contain power values for part of a power spectrum)
    """
    # # Preload variable to indicate the length of the total spectrum (all intervals together)
    length_of_spectrum = 0

    # # Add interval lengths to calculate length of total spectrum
    for i in range(0, len(y_intervals)):
        length_of_spectrum += len(y_intervals[i])

    # # Preload list for calculating signal/noise ratio
    y_over_noise = []

    # # Loop to calculate signal/noise ratio
    for i in range(0, len(y_intervals)):
        if median[i] == 0: # # Check to avoid error if median is zero (as this is not necessarily a problem).
            median[i] = 0.000001
            print('Dividing by zero')

        # # Calculate signal/noise for current signal value
        y_temp = [y / median[i] for y in y_intervals[i]]

        # # Add signal/noise temporary list to end of result list y_over_noise
        y_over_noise = y_over_noise + y_temp
    return y_over_noise


def good_peak_finder(peaks_over_noise, p_peaks, freq_peaks, freq_muHz_peaks, intervals, pn_limit):
    # # Function to find good peaks by comparing peak height with noise, and requiring a signal/noise ratio over
    # 4-6 (changed manually in function). Can receive multiple intervals (x-axis start and end) and corresponding data.
    """
    Function that, given peaks found by peak-finding, and knowing their matching signal/noise ratios, can sort them
    based on a pre-determined signal/noise lower limit (it checks whether each peak is above this limit, and rearranges
    in new lists). The input variables are structured in this sense, because it expects to receive variables from
    func_over_noise and median_noise.


    Variables:
        peaks_over_noise: List of signal/noise values for all the peaks found by peak-finding (all intervals).
                          Structured based on the variable y_over_noise from func_over_noise.

        p_peaks: List of lists, with peak power for each interval. Similar in structure as y_intervals from
                 median_noise, just with less elements in each interval (only peaks).

        freq_peaks: List of lists, same as p_peaks, but with frequency value for each peak instead of peak power.

        freq_muHz_peaks: List of lists, same as freq_peaks, but values are in microHz.

        intervals: List of lists, with each list including start and finish integer for the corresponding interval.
                   Similar in structure to variable intervals from median_noise

        pn_limit: pre-specified power/noise lower limit (f.ex. 6). Float or integer


    Returns:
        [goodstuff, totalstuff]
            goodstuff: List of lists, [good_peaks, good_freq_muHz, good_freq, good_peaks_over_noise], that includes info
                       on the peaks above the power/noise limit only.
                       good_peaks: list of floats, with power of all the peaks above p/n limit.
                       good_freq_muHz & good_freq: list of floats, with frequency (in microHz and Hz, respectively) of
                                                   the peaks found above the p/n limit
                       good_peaks_over_noise: list of floats, with power/noise level of all the peaks above p/n limit.

            totalstuff: List of lists, [p_peaks_tot, freq_muHz_peaks_tot, freq_peaks_tot], that includes info of all the
                        peaks, both above and below power/noise limit. Similar to goodstuff in structure, but missing
                        a list of power/noise levels.

    """
    # # Preload lists for full spectrum (not specific intervals)
    p_peaks_tot = []
    freq_peaks_tot = []
    freq_muHz_peaks_tot = []

    # # Loop to add all intervals' data to the full spectrum lists
    for i in range(0, len(intervals)):
        p_peaks_tot = p_peaks_tot + p_peaks[i]
        freq_peaks_tot = freq_peaks_tot + freq_peaks[i]
        freq_muHz_peaks_tot = freq_muHz_peaks_tot + freq_muHz_peaks[i]

    # # Preload lists that only includes peaks with a signal/noise ratio above the indicated
    good_peaks = []
    good_peaks_over_noise = []
    good_freq = []
    good_freq_muHz = []

    # # Loop to find peaks to include in the good peaks lists
    for i in range(0, len(peaks_over_noise)):
        if peaks_over_noise[i] > pn_limit:  # in this check, the signal/noise limit is stated
            # # Append to pre-loaded lists
            good_peaks.append(p_peaks_tot[i])
            good_peaks_over_noise.append(peaks_over_noise[i])
            good_freq.append(freq_peaks_tot[i])
            good_freq_muHz.append(freq_muHz_peaks_tot[i])

    # # Group in return lists
    goodstuff = [good_peaks, good_freq_muHz, good_freq, good_peaks_over_noise]
    totalstuff = [p_peaks_tot, freq_muHz_peaks_tot, freq_peaks_tot]
    return [goodstuff, totalstuff]


def autocorrellator(freq, spectrum):  # Depreciated, should only be used for testing autocorr function
    """
    Function to calculate autocorrelation of power-spectrum using the mathematical definition of autocorrelation.
    Very slow method. Not used in analysis, depreciated and replaced by the function autocorr.

    Returns:
        f_change: List of floats, giving displacement frequency for each autocorrelation element.
        autocorrelation: List of floats, giving autocorrelation value for each displacement frequency.

    Variables:
        freq: List of floats, with frequency values for the whole power-spectrum. Also used to define the step size for
              this autocorrelation function.
        spectrum: List of floats, each element gives power of the spectrum, matching the frequency values in freq.
    """

    # Create list with frequency used for reference, and frequencies to be displaced for correlation
    reference_frequency = freq
    correlation_frequency = [x for x in freq]

    # Introduce padding-weights around ends
    weights = [np.sin(np.pi*x/len(spectrum)) for x in range(0, len(spectrum))]
    spectrum_weighted = [spectrum[x] * weights[x] for x in range(0, len(spectrum))]

    # Find step-size for correlation
    d_f = reference_frequency[1] - reference_frequency[0]

    n = len(reference_frequency)

    # Pre-load lists of displacement frequency and autocorrelation
    f_change = [0] * int(n/5)
    autocorrelation = [0] * int(n/5)

    # Loop to calculate autocorrelation
    for i in range(0, int(n/5)):

        # Loop to displace correlation frequency by d_f for each element, and loop around if it reaches the end of
        # reference_frequency
        for j in range(0, len(correlation_frequency)):
            correlation_frequency[j] = correlation_frequency[j] - d_f
            if correlation_frequency[j] < reference_frequency[0]:
                correlation_frequency[j] = correlation_frequency[j] + reference_frequency[-1]

        # Add current displacement frequency to f_change list
        if i > 0:
            f_change[i] = f_change[i-1] + d_f
        elif i == 0:
            f_change[i] = d_f
        else:
            print('Error')

        # Preload temporary list of products between displaced and undisplaced function values, to be summed together
        multiplication_list = [0]*len(reference_frequency)

        # Loop over elements of function
        for k in range(0, len(reference_frequency)):

            # Find the index, where the correlation frequency matches reference_frequency[0] (as such, find the index
            # difference between the two for the current step)
            if k == 0:
                # Float rounding means the correlation_frequency might be slightly different from reference_frequency
                try:
                    correlation_index = correlation_frequency.index(reference_frequency[k])
                except ValueError:
                    # Slow method, find nearest match. Do this only if rounding error made previous attempt
                    # unsuccessful
                    abs_diff = [abs(correlation_frequency[x] - reference_frequency[k]) for x in range(0, n)]
                    correlation_index = np.argmin(abs_diff)

            # Loop around to beginning if end is reached
            elif correlation_index == len(reference_frequency)-1:
                correlation_index = 0

            # Happens if correlation_index is already calibrated to match reference_frequency. As such, just advance one
            # step (same as reference_frequency)
            else:
                correlation_index += 1

            # Calculate the correlation between the two functions at the current index
            try:
                multiplication_list[k] = spectrum_weighted[k] * spectrum_weighted[correlation_index]
            except IndexError:
                print(len(multiplication_list))
                print(k)
                print(len(reference_frequency))
                print(correlation_index)

        # Sum together multiplication_index to find the autocorrelation for current displacement frequency
        autocorrelation[i] = np.sum(multiplication_list)

    # Plot correlation to assess
    plt.plot(f_change, autocorrelation)
    plt.show()

    return [f_change, autocorrelation]


def autocorr(x):
    """
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.correlate.html

    Function to calculate autocorrelation, created with inspiration from http://stackoverflow.com/q/14297012/190597.
    Uses the numpy function "correlate", that calculates the correlation as a convolution.
    The frequency displacement is not needed to calculate the autocorrelation, but can be deduced by knowing the
    difference in frequency between each element of x.

    Returns:
        result: a list of normed autocorrelation values, one for each displacement frequency

    Variables:
        x: List of float, signal data usually, in this context it would be the power spectrum values
    """

    # # Convert list to numpy array
    data = np.asarray(x)

    # # Calculate sum squared, to norm autocorrelation
    norm = np.sum(data**2)

    # # Calculate correlation from negative to positive frequency change
    res = np.correlate(data, data, mode='full')

    # # Save only the positive part (last half)
    res = res[res.size//2:]

    # # Norm data by variance
    result = res/norm
    return result


def echelle2(x, modulus, **kwargs):
    # # Allows weights to be used to create a heatmap. Simply give the function the optional argument weights=some_list,
    # where some_list is a list of len(x), of floats or integers that is to be used as a weight when creating the
    # heatmap. A possible use is peak height.

    """
    A function used for plotting echelle diagrams, and using them to calculate small frequency separations.
    Called echelle2 only because there once existed a function called echelle, that had less optional inputs.

    Returns:
        Nothing, or
        points (list of x and y plot coordinates to be used to calculate small frequency separation), or
        [fig, ax] (figure and axes handle, only if called in animation mode)

    Variables:
        x: List of floats, intended to be frequency values.
        modulus: Float, value used to calculate modulo of x (x mod modulus).

    Optional variables:
        heatmap: Boolean, determines whether a heatmap is plotted, otherwise, a normal point-plot occurs.
        weights: List of floats, used for heatmap colouring to determine height of peaks.
        degree_coords: List of 3 lists, with xy-coordinates used to plot angular order degree numbers.
        plot_points: List of lists, xy-coordinates to plot previously selected peaks for small freq separation calc.
        animation_mode: Boolean, determines whether figure and axes handle should be returned for looping
                        (Used to animate the plot for different modulus values)
        figure: List [fig, ax], with figure and axes handle
        xlim: List, with x-axis limits
        ylim: List, with y-axis limits
        xlab: String, label for x-axis
        ylab: String, label for y-axis
        mark_mode: Boolean, determines whether graphical input is used to mark placement of angular degree
                   numbers (for later plots), or to mark points used to calculate small frequency separation.
        block: Boolean, used to determine whether plot should block advancement of code.
        nbins: Integer, determines number of heatmap bins. Default 64
    """
    # # Get boolean to decide whether heatmap should be created
    heatmap = kwargs.get('heatmap', False)

    # # Get optional weights for heatmap
    weights = kwargs.get('weights', [1] * len(x))

    # # Get optional coordinates to plot angular degree numbers
    degree_coords = kwargs.get('degree_coords', None)

    # # Get optional coordinates to plot points
    plot_points = kwargs.get('plot_points', None)

    # # Get boolean to decide if the function should be run in animation mode
    animation_mode = kwargs.get('animation_mode', False)

    # # Get figure (used in animation mode)
    [fig, ax] = kwargs.get('figure', plt.subplots())

    # # Get x and y axis limit
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)

    # # Get optional labels for plot axis
    ylab = kwargs.get('ylabel', None)
    xlab = kwargs.get('xlabel', None)

    # # Get optional marking mode
    mark_mode = kwargs.get('mark_mode', False)

    # # Get optional block mode
    block = kwargs.get('block', True)

    # # Calculate x mod modulus
    xmod = [math.fmod(i, modulus) for i in x]

    # # Plot Echelle diagram
    if heatmap is False:

        # Check and run animation mode

        # Make plot
        ax.plot(xmod, x, 'r*')

        # Make title
        ax.set_title('Echelle Diagram')

        # Check for labels
        if xlab is not None:
            ax.set_xlabel(xlab)
        if ylab is not None:
            ax.set_ylabel(ylab)

        # Check for limits
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        # Check for animation mode and show plot
        if animation_mode is True:
            pass
        else:
            plt.show()

    # # Optional, plot heatmap Echelle diagram
    if heatmap is True:

        # Get number of bins
        nbins = kwargs.get('number_of_bins', 64)

        # Construct 2D histogram from data using the 'plasma' colormap
        ax.hist2d(xmod, x, bins=nbins, weights=weights, normed=False, cmap='plasma')

        # Make title
        ax.set_title('Heatmap Echelle diagram')

        # Check for labels
        if xlab is not None:
            ax.set_xlabel(xlab)
        if ylab is not None:
            ax.set_ylabel(ylab)

        # Check for limits
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        # Check for degree labels and plot them
        if degree_coords is not None:
            for i in range(0, len(degree_coords)):
                current_degree = degree_coords[i]
                ax.text(current_degree[0], current_degree[1], str(i), color='w')

        # Check for points and plot them
        if plot_points is not None:
            for i in range(0, len(plot_points)):
                current_point = plot_points[i]
                ax.plot(current_point[0], current_point[1], 'rX', markersize=4)

        # Check for animation mode and show plot
        if animation_mode is True:
            pass
        elif mark_mode is True:
            points = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
        else:
            plt.show(block=block)

    if animation_mode is True:
        return [fig, ax]
    elif mark_mode is True:
        return points


def gauss_smooth(y, window_size=3, sigma=2):
    """
    source: https://github.com/USGS-Astrogeology/PyHAT/blob/master/libpysat/spectral/smoothing.py
    Apply a gaussian filter to smooth the input vector

    Parameters
    ==========
    y :  array
         The input array

    window_size : int
                  An odd integer describing the size of the filter.

    sigma : float
            The number of standard deviation
    """
    filt = signal.windows.gaussian(window_size, sigma, sym=False)
    print(window_size)
    return Series(signal.convolve(y, filt, mode='same'), index=y.index)


def hann_smooth(y, window_size=3, scale=True):
    """
    source: https://github.com/USGS-Astrogeology/PyHAT/blob/master/libpysat/spectral/smoothing.py
    Apply a gaussian filter to smooth the input vector

    Parameters
    ==========
    y :  array
         The input array

    window_size : int
                  An odd integer describing the size of the filter.

    sigma : float
            The number of standard deviation
    """
    filt = signal.windows.hann(window_size, sym=False)
    print(window_size)
    res = Series(signal.convolve(y, filt, mode='same'), index=y.index)
    if scale is True:
        return res/window_size
    else:
        return res


def smoother(freq, p, scale=True, scale_man=True, window_size=False):
    """
    Function that prepares signal data to be smoothed, and smoothes data with a Hann window function by referencing
    hann_smooth.

    Returns:
        p_smooth: List of floats, smoothed spectrum values

    Variables:
        freq: List of floats, frequency values for the spectrum
        p: List of floats, y-values for spectrum (can be power or power/noise or similar)

    Optional variables:
        scale: Boolean, determines whether or not the smoothed spectrum should be scaled by dividing with the window
               size. Default is True.
        scale_man: Boolean, determines whether or not the smoothed spectrum should be scaled to have the same max value
                   as the original spectrum. Useful for plotting. Default is True.
        window_size: Float or Integer. Makes it possible to set the window size when calling the
                     function. Otherwise, a prompt will be issued.
    """

    # # Plot and select smoothing window size (in frequency size)
    if window_size is False:
        plt.plot(freq, p)
        plt.show()
        wsize_freq = float(input('What frequency size should the window have?'))
    else:
        wsize_freq = window_size

    # Frequency per index
    freq_ind = np.mean(np.diff(freq))

    # Window size in index, and set to nearest odd integer (rounding up if equal)
    wsize = wsize_freq / freq_ind
    if np.fmod(int(np.ceil(wsize)), 2) == 0:
        if np.floor(wsize) == np.ceil(wsize):
            wsize = int(np.floor(wsize)) + 1
        else:
            wsize = int(np.floor(wsize))
    elif np.fmod(int(np.floor(wsize)), 2) == 0:
        if np.ceil(wsize) == np.floor(wsize):
            wsize = int(np.ceil(wsize)) + 1
        else:
            wsize = int(np.ceil(wsize))
    else:
        print('Error in finding odd wsize')
        wsize = int(wsize)

    print('wsize', wsize)

    # # Smooth
    p_array = Series(p)
    p_smooth = hann_smooth(p_array, window_size=wsize, scale=scale).tolist()
    maxp = np.max(p)
    maxs = np.max(p_smooth)
    if scale_man is True:
        p_smooth[:] = [maxp * x / maxs for x in p_smooth]
    plt.plot(freq, p, 'C6--', linewidth=0.5)
    plt.plot(freq, p_smooth, linewidth=1)
    plt.show()

    return p_smooth
