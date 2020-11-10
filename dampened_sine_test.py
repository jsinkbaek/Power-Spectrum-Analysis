import numpy as np
import matplotlib.pyplot as plt
import create_pspectrum
import clean_procedure


def gen_noisy_sinusoid(timestamps, freq, phase):
    noise = np.random.normal(0, 0.05, timestamps.size)
    sinusoid = np.sin(2*np.pi*timestamps*freq + phase) + noise
    return sinusoid


def gen_noisy_damped_sinusoid(timestamps, freq, phase, decay):
    noise = np.random.normal(0, 0.015, timestamps.size)
    sinusoid = np.exp(-decay*timestamps) * np.sin(2*np.pi*timestamps*freq + phase) + noise
    return sinusoid


def gen_timestamps(t0, t1, spacing):
    """
    Generates semi-equispaced timestamps, with random points taken out
    :param t0: starting time
    :param t1: end time
    :param spacing: equispaced interval
    :return: equispaced timestamps, but with random points removed
    """
    steps = int((t1 - t0)/spacing)
    timestamps = np.linspace(t0, t1, steps)
    rm_amount = np.random.randint(timestamps.size / 10, timestamps.size / 5)
    deselect = np.unique(np.random.randint(0, timestamps.size, size=rm_amount))
    timestamps = np.delete(timestamps, deselect)
    return timestamps


signal_freq = 20
freq_centre = [21]
t = gen_timestamps(0, 1, 0.001)

signal1 = gen_noisy_sinusoid(t, freq=signal_freq, phase=0)
plt.figure()
plt.plot(t, signal1, 'r')
plt.title('Undampened sine wave, freq 0.05, phase 0')
plt.show(block=False)
frequency, power, alpha, beta = create_pspectrum.cuda(signal1, t, freq_centre=freq_centre, half_width=signal_freq,
                                                      resolution=0.0001, chunk_size=1000, silent=True)[0]
plt.figure()
plt.plot(frequency, power)
plt.title('Undampened sine wave power')
plt.show(block=False)

signal2 = gen_noisy_sinusoid(t, freq=signal_freq, phase=0.34)
plt.figure()
plt.plot(t, signal1, 'r', t, signal2, 'b', t, signal1+signal2, 'k--')
plt.title('2 undampened sine waves, freq 0.05')
plt.legend(['phase 0', 'phase 0.34', 'sum'])
plt.show(block=False)
frequency, power, alpha, beta = create_pspectrum.cuda(signal1+signal2, t, freq_centre=freq_centre, half_width=signal_freq,
                                                      resolution=0.0001, chunk_size=1000, silent=True)[0]
plt.figure()
plt.plot(frequency, power)
plt.title('Sum of two undampened sine waves')
plt.show(block=False)

ts1 = gen_timestamps(0, 0.5, 0.001)
ts2 = gen_timestamps(0.5, 1, 0.001)
signal3 = gen_noisy_sinusoid(ts1, freq=signal_freq, phase=0)
signal4 = gen_noisy_sinusoid(ts2, freq=signal_freq, phase=0.34)
ts12 = np.append(ts1, ts2)
signal34 = np.append(signal3, signal4)
plt.figure()
plt.plot(ts1, signal3, 'r', ts2, signal4, 'b')
plt.title('two signals freq 0.05, different phase, stops abruptly')
plt.legend(['phase 0', 'phase 0.34'])
plt.show(block=False)
frequency, power, alpha, beta = create_pspectrum.cuda(signal34, ts12, freq_centre=freq_centre, half_width=signal_freq,
                                                      resolution=0.0001, chunk_size=1000, silent=True)[0]
plt.figure()
plt.plot(frequency, power)
plt.title('Two abruptly ending sine waves')
plt.show(block=False)

ts3 = gen_timestamps(0, 1, 0.001)
signal5 = gen_noisy_damped_sinusoid(ts3, freq=signal_freq, phase=0, decay=8)
plt.figure()
plt.plot(ts3, signal5)
plt.title('One dampened sine wave')
plt.show(block=False)
frequency, power, alpha, beta = create_pspectrum.cuda(signal5, ts3, freq_centre=freq_centre, half_width=signal_freq,
                                                      resolution=0.0001, chunk_size=1000, silent=True)[0]
plt.figure()
plt.plot(frequency, power)
plt.title('One dampened sine wave')
plt.show(block=False)
signal6 = gen_noisy_damped_sinusoid(ts3, freq=signal_freq, phase=0.34, decay=4)
ts4 = 1 - ts3       # reverse timestamps for mirrored decaying sinusoids
signal56 = signal5 + signal6[::-1]
plt.figure()
plt.plot(ts3, signal5, 'r', ts4, signal6, 'b', ts3, signal56, 'k--')
plt.title('two dampened waves (one mirrored)')
plt.show(block=False)
frequency, power, alpha, beta = create_pspectrum.cuda(signal56, ts3, freq_centre=freq_centre, half_width=signal_freq,
                                                      resolution=0.0001, chunk_size=1000, silent=True)[0]
plt.figure()
plt.plot(frequency, power)
plt.title('Two dampened sine waves')
plt.show(block=True)
clean_procedure.cuda(t, signal1, 4, freq_centre=freq_centre[0], half_width=signal_freq, resolution=0.0001,
                     silent=False, mph=0.05)
clean_procedure.cuda(ts3, signal5, 20, freq_centre=freq_centre[0], half_width=signal_freq, resolution=0.0001,
                     silent=False, mph=0.0001)
clean_procedure.cuda(ts3, signal56, 20, freq_centre=freq_centre[0], half_width=signal_freq, resolution=0.0001,
                     silent=False, mph=0.0001)
