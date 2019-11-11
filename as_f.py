# # # # # # #
# This is a function "repository", that stores some functions used for the asterseismology-scripts.
# # # # # # #

# # Import Statements for functions
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import signal
from pandas import Series


def autocorr(y, number_of_steps=None, x=None, x_tot=None):
    """
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.correlate.html

    Function to calculate autocorrelation, created with inspiration from http://stackoverflow.com/q/14297012/190597.
    Uses the numpy function "correlate", that calculates the correlation as a convolution.
    The frequency displacement is not needed to calculate the autocorrelation, but can be deduced by knowing the
    difference in frequency between each element of x.

    Returns:
        result: a list of normalized autocorrelation values, one for each displacement frequency

    Variables:
        y: List of float, signal data usually, in this context it would be the power spectrum values
    """
    # # Convert lists to numpy array
    y = np.asarray(y)
    x = np.asarray(x)
    x_tot = np.asarray(x_tot)

    # # Create autocorrelation from smaller data set, while maintaining same steps as defined by number_of_steps
    if number_of_steps is not None:
        if len(y) != number_of_steps and len(y) == len(x) and len(x_tot) == number_of_steps:
            data = np.zeros(number_of_steps)
            for i in range(0, len(x)):
                indx = np.argmin(np.abs(x_tot - x[i]))
                data[indx], x_tot[indx] = y[i], x[i]
        else:
            data = y
            print('Warning, number_of_steps specified, but problems with y, x, or x_tot.')
    else:
        data = y

    # # Calculate sum squared, to norm autocorrelation
    norm = np.sum(data**2)

    # # Calculate correlation from negative to positive frequency change
    res = np.correlate(data, data, mode='full')

    # # Save only the positive part (last half)
    res = res[res.size//2:]

    # # Norm data by sum
    result = res/norm
    if x_tot is not None:
        return result, x_tot
    else:
        return result


def echelle(x, modulus, **kwargs):
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
    Apply a hann filter to smooth the input vector

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
