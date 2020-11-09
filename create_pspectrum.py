import math
import numpy as np
import time as tm
from detect_peaks import detect_peaks
import matplotlib.pyplot as plt


def numpy(y, t, freq_centre, half_width, resolution, chunk_size=50, dtype=np.double):
    """
    Calculate powerspectrum using matrix multiplication with numpy
    """
    pi = math.pi
    t1 = tm.time()
    # # Prepare input for use

    # Convert input to numpy array (if not already)
    y = np.asarray(y)
    t = np.asarray(t)

    # Change to angular frequencies (assumes input is in cyclic frequencies)
    freq_centre[:] = [x * (2 * pi) for x in freq_centre]
    half_width = half_width * 2 * pi
    resolution = resolution * 2 * pi

    # Data mean subtraction, to reduce potentially constant elements
    y_mean = np.mean(y, dtype=dtype)
    y = y - y_mean
    t = t - np.min(t)  # To get first element in t to be 0

    # # Prepare for for-loop

    # Preload results list
    results = [None] * len(freq_centre)

    # Set amount of steps (might be subject to rounding error)
    step_amnt = int((2 * half_width) / resolution)

    # # # # Run for loop for all frequencies indicated by length of freq_centre # # # #
    for k in range(0, len(freq_centre)):
        # Show numpy config to check which BLAS is used
        # np.__config__.show()

        # Get current frequency centre
        freq_centre_current = freq_centre[k]

        # Create frequency steps
        freq = np.linspace(freq_centre_current - half_width, freq_centre_current + half_width, step_amnt, dtype=dtype)

        # Reshape freq and t in order to do matrix multiplication
        freq = np.reshape(freq, (len(freq), 1))
        lent = len(t)

        # # Zhu Li, do the thing

        # Calculate sine and cosine function part values
        print(freq.shape)
        print(t.shape)

        sin = np.zeros(step_amnt)
        cos = np.zeros(step_amnt)
        sin2 = np.zeros(step_amnt)
        cos2 = np.zeros(step_amnt)
        sincos = np.zeros(step_amnt)

        freq = np.ascontiguousarray(freq)
        t = np.ascontiguousarray(t)

        # Recurrence sine and cosine difference product
        s_diff = np.sin(resolution * t, dtype=dtype)
        c_diff = np.cos(resolution * t, dtype=dtype)

        # # # # # Calculation matrices # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # use [c0, s0][c_diff, s_diff; -s_diff, c_diff]  (inverted in calc for .T)         #
        # Recurrence based on s_m = c_(m-1)*sin(deltaf * t) + s_(m-1)*cos(deltaf * t) and  #
        # c_m = c_(m-1)*cos(deltaf * t) - s_(m-1)*sin(deltaf * t) from T. Ponman 1981      #
        # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        calc_base = np.array([[c_diff, -s_diff], [s_diff, c_diff]]).T
        print('calc_base.shape', calc_base.shape)
        calc_mat = np.zeros((chunk_size, lent, 2, 2))
        calc_mat[0, :, :, :] = calc_base

        # Generates the necessary matrices (multiplied) for each frequency point based on the original point
        for i in range(1, chunk_size):
            calc_mat[i, :, :, :] = np.matmul(calc_mat[i - 1, :, :, :], calc_base)

        # Convert large matrix arrays to contigous arrays
        calc_mat = np.ascontiguousarray(calc_mat)

        # Calculates sine and cosine values
        for i in range(0, step_amnt, chunk_size):
            end = i + chunk_size
            if end > step_amnt:
                dif = end - step_amnt
                calc_mat = np.delete(calc_mat, range(chunk_size-dif, chunk_size), 0)
                end = step_amnt
                chunk_size = end - i

            print('Current chunk ', i, ':', end, ' of ', step_amnt)

            # Original point calculation
            s0 = np.sin(freq[i] * t, dtype=dtype)
            c0 = np.cos(freq[i] * t, dtype=dtype)

            # Sine/cosine vector initialization (for matmul calculation, see c0, s0 before loop)
            trig_vec = np.zeros((1, lent, 1, 2))
            trig_vec[0, :, 0, 0] = c0
            trig_vec[0, :, 0, 1] = s0
            trig_vec = np.repeat(trig_vec, chunk_size, axis=0)

            # Matrix calculations)
            matrix_result = np.matmul(trig_vec, calc_mat)
            sin_temp = matrix_result[:, :, 0, 1]
            cos_temp = matrix_result[:, :, 0, 0]

            # Sum and save results
            sin[i:end] = np.sum(y * sin_temp, 1)
            cos[i:end] = np.sum(y * cos_temp, 1)
            sin2[i:end] = np.sum(sin_temp ** 2, 1)
            cos2[i:end] = np.sum(cos_temp ** 2, 1)
            sincos[i:end] = np.sum(sin_temp * cos_temp, 1)

        # # Calculate alpha and beta components of spectrum, and from them, power of spectrum
        alpha = (sin * cos2 - cos * sincos) / (sin2 * cos2 - sincos**2)
        beta = (cos * sin2 - sin * sincos) / (sin2 * cos2 - sincos**2)

        power = alpha**2 + beta**2

        # # Last loop steps

        # Convert frequency back to cyclic units
        freq = freq / (2 * pi)

        # Save data in results
        results[k] = [freq, power, alpha, beta]
    t2 = tm.time()
    print('Total time elapsed: ', t2-t1, ' seconds')

    return results


def cuda(y, t, freq_centre, half_width, resolution, chunk_size=100, dtype=None):
    """
    Calculate powerspectrum using matrix multiplication with cupy (CUDA numpy)
    """
    # # Uses cupy (cuda numpy) instead of numpy # #
    import cupy as cu
    if dtype is None:
        dtype = cu.double

    pi = math.pi
    t1 = tm.time()

    # # # # Prepare input for use # # # #
    # Convert input to cupy array (if not already)
    y = cu.asarray(y)
    t = cu.asarray(t)

    # Change to angular frequencies (assumes input is in cyclic frequencies)
    freq_centre[:] = [x * (2 * pi) for x in freq_centre]
    half_width = half_width * 2 * pi
    resolution = resolution * 2 * pi

    # Data mean subtraction, to reduce potentially constant elements
    y_mean = cu.mean(y, dtype=dtype)
    y = y - y_mean
    t = t - cu.min(t)  # To get first element in t to be 0

    # # # # Prepare for for-loop # # # #
    # Preload results list
    results = [None] * len(freq_centre)

    # Set amount of steps (might be subject to rounding error)
    step_amnt = int((2 * half_width) / resolution)

    # # # # Run for loop for all frequencies indicated by the list freq_centre # # # #
    for k in range(0, len(freq_centre)):
        freq_centre_current = freq_centre[k]
        # Create frequency steps
        freq = cu.linspace(freq_centre_current - half_width, freq_centre_current + half_width, step_amnt, dtype=dtype)
        # Reshape freq in order to do matrix multiplication
        freq = cu.reshape(freq, (len(freq), 1))
        lent = len(t)

        # # Prepare to calculate sine and cosine function parts of power spectrum estimate # #
        sin = np.zeros(step_amnt)
        cos = np.zeros(step_amnt)
        sin2 = np.zeros(step_amnt)
        cos2 = np.zeros(step_amnt)
        sincos = np.zeros(step_amnt)

        freq = cu.ascontiguousarray(freq)
        t = cu.ascontiguousarray(t)
        # Recurrence sine and cosine difference product
        s_diff = cu.sin(resolution * t, dtype=dtype)
        c_diff = cu.cos(resolution * t, dtype=dtype)

        # # # # # Calculation matrices # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # use [c0, s0][c_diff, s_diff; -s_diff, c_diff]  (inverted in calc for .T)         #
        # Recurrence based on s_m = c_(m-1)*sin(deltaf * t) + s_(m-1)*cos(deltaf * t) and  #
        # c_m = c_(m-1)*cos(deltaf * t) - s_(m-1)*sin(deltaf * t) from T. Ponman 1981      #
        # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        # # Prepare linear operation matrix with initial operation # #
        calc_base = cu.array([[c_diff, -s_diff], [s_diff, c_diff]]).T
        calc_mat = cu.zeros((chunk_size, lent, 2, 2))
        calc_mat[0, :, :, :] = calc_base
        # # Generate linear operation matrices for recurrence multiplication # #
        for i in range(1, chunk_size):
            calc_mat[i, :, :, :] = cu.matmul(calc_mat[i - 1, :, :, :], calc_base)

        # Convert large matrix arrays to contiguous arrays
        calc_mat = cu.ascontiguousarray(calc_mat)

        # # Calculate sine and cosine function parts of power spectrum estimation # #
        for i in range(0, step_amnt, chunk_size):
            end = i + chunk_size
            if end > step_amnt:
                dif = end - step_amnt
                calc_mat = cu.asnumpy(calc_mat)
                calc_mat = np.delete(calc_mat, range(chunk_size - dif, chunk_size), 0)
                calc_mat = cu.array(calc_mat)
                calc_mat = cu.ascontiguousarray(calc_mat)
                end = step_amnt
                chunk_size = end - i

            print('Current chunk ', i, ':', end, ' of ', step_amnt)
            # Original point calculation
            s0 = cu.sin(freq[i] * t, dtype=dtype)
            c0 = cu.cos(freq[i] * t, dtype=dtype)
            # Sine/cosine vector initialization (for matmul calculation, see c0, s0 before loop)
            trig_vec = cu.zeros((1, lent, 1, 2))
            trig_vec[0, :, 0, 0] = c0
            trig_vec[0, :, 0, 1] = s0
            trig_vec = cu.repeat(trig_vec, chunk_size, axis=0)
            # Matrix calculations
            matrix_result = cu.matmul(trig_vec, calc_mat)
            sin_temp = matrix_result[:, :, 0, 1]
            cos_temp = matrix_result[:, :, 0, 0]

            # # Sum and save results # #
            sin[i:end] = cu.sum(y * sin_temp, 1).get()
            cos[i:end] = cu.sum(y * cos_temp, 1).get()
            sin2[i:end] = cu.sum(sin_temp ** 2, 1).get()
            cos2[i:end] = cu.sum(cos_temp ** 2, 1).get()
            sincos[i:end] = cu.sum(sin_temp * cos_temp, 1).get()

        # # # # Calculate alpha and beta components of spectrum, and from them, power of spectrum # # # #
        alpha = (sin * cos2 - cos * sincos) / (sin2 * cos2 - sincos**2)
        beta = (cos * sin2 - sin * sincos) / (sin2 * cos2 - sincos**2)
        power = alpha**2 + beta**2

        # # # # Last loop steps # # # #
        # Convert frequency back to cyclic units
        freq = freq.get() / (2 * pi)
        # Save data in results
        results[k] = [freq, power, alpha, beta]
    t2 = tm.time()
    print('Total time elapsed create_pspectrum_cuda: ', t2-t1, ' seconds')
    return results


def nufftpy(y, t, half_width, resolution):
    import nufftpy as nfpy
    steps = int((2*half_width)/resolution)
    freq = nfpy.nufftfreqs(steps, df=resolution)
    freq = freq[len(freq)//2:-1]
    harmonic_content = nfpy.nufft1(t, y, steps, df=(resolution*2*math.pi))
    harmonic_content = harmonic_content[len(harmonic_content)//2:-1]

    return [freq, harmonic_content.real**2+harmonic_content.imag**2, harmonic_content.real, harmonic_content.imag]
