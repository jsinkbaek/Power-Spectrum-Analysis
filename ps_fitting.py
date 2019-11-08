import matplotlib.pyplot as plt
import ps_f
from scipy.optimize import curve_fit
from detect_peaks import detect_peaks
import lmfit as lm
import numpy as np

print('Which spectrum to plot? Insert datafile, excluding .dat')
inpt1 = input()

freq = []
p = []
data = ps_f.preader(inpt1 + '.dat')
for row in data:
    freq.append(row[0])
    p.append(row[1])

muHz = 0.000001
freq_muHz = [x / muHz for x in freq]


def fitting_double_interval(freq, p):
    plt.plot(p)
    plt.title('Select end of first interval, start of second interval, Lorentz fit')
    coords = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    coord1 = 0
    coord2 = coords[0]
    coord2 = int(coord2[0])
    coord3 = coords[1]
    coord3 = int(coord3[0])
    coord4 = len(p)
    plt.close()

    p_fit = p[coord1:coord2] + p[coord3:coord4]
    freq_fit = freq[coord1:coord2] + freq[coord3:coord4]

    from lmfit.models import ConstantModel, LorentzianModel
    const_mod = ConstantModel()
    lorz_mod = LorentzianModel()

    pars = const_mod.make_params()
    pars += lorz_mod.make_params()

    mod = const_mod + lorz_mod
    output = mod.fit(p_fit, pars, x=freq_fit)
    dicc = output.best_values  # dictionary of fit values
    print(output.fit_report())
    print(output.best_values)
    c, amplitude, sigma, center = dicc['c'], dicc['amplitude'], dicc['sigma'], dicc['center']
    param_res = [c, amplitude, sigma, center]

    plt.plot(freq_fit, p_fit, 'r*')
    fit = [ps_f.lorentz2(x, amplitude, sigma, center, c) for x in freq]
    plt.plot(freq_fit, output.best_fit)
    plt.xlim(freq_fit[0], freq_fit[-1])
    plt.show()

    p_res1 = [i - j for i, j in zip(p, fit)]
    for k in range(0, len(p_res1)):
        if p_res1[k] <= 0:
            p_res1[k] = 0
    plt.plot(freq, p_res1)
    plt.show()
    return[p_res1, param_res]


def fitting_peaks_double_interval(freq, p, MPH):
    plt.plot(p)
    x_fit_values = list(range(0, len(p)))
    plt.title('Select end of first fit interval, and start of second.')
    coords = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)
    coord1 = 0
    coord2 = coords[0]
    coord2 = int(coord2[0])
    coord3 = coords[1]
    coord3 = int(coord3[0])
    coord4 = len(p)
    plt.close()

    p_temp = p[coord1:coord2] + p[coord3:coord4]
    x_fit_values_temp = x_fit_values[coord1:coord2] + x_fit_values[coord3:coord4]
    freq_temp = freq[coord1:coord2] + freq[coord3:coord4]

    indices_peaks = detect_peaks(p_temp, mph=MPH, show=True)
    peaks = [0] * len(indices_peaks)
    freq_peaks = [0] * len(indices_peaks)
    x_peaks = [0] * len(indices_peaks)
    for l in range(0, len(indices_peaks)):
        peaks[l] = p_temp[indices_peaks[l]]
        freq_peaks[l] = freq_temp[indices_peaks[l]]
        #x_peaks[l] = x_fit_values_temp[indices_peaks[l]]
    x_peaks = freq_peaks
    from lmfit.models import ConstantModel, LorentzianModel
    from lmfit import Model
    const_mod = ConstantModel()
    lorz_mod = LorentzianModel()


    pars = const_mod.make_params()
    pars += lorz_mod.make_params()

    mod = const_mod + lorz_mod
    print(pars)
    #pars['center'].set(value=-5000, max=-1)
    #pars['amplitude'].set(value=2000000, min=1)
    #pars['c'].set(value=10, max=25)
    #pars['sigma'].set(value=4000, min=100)
    #pars['center'].set(value=-0.005, max=-0.000001, min=-1)
    #pars['amplitude'].set(value=200, min=200)
    pars['c'].set(value=0, max=0.000000001, min=0)
    #pars['sigma'].set(value=4000, min=100)
    output = mod.fit(peaks, pars, x=x_peaks)
    print(output.fit_report())
    print(output.best_values)
    #print(output.best_fit)
    dicc = output.best_values  # dictionary of fit values
#    for key,val in dicc.items():  # Should convert into variables, does not
#        exec(key + '=val')
    c, amplitude, sigma, center = dicc['c'], dicc['amplitude'], dicc['sigma'], dicc['center']
    param_res = [c, amplitude, sigma, center]
    fit_data = [ps_f.lorentz2(x, amplitude, sigma, center, c) for x in x_peaks]
    plt.plot(x_peaks, peaks, 'r*')
    plt.plot(x_peaks, output.best_fit, 'b-')
    plt.show()

    p_res2 = [i - j for i, j in zip(p, fit_data)]
    for h in range(0, len(p_res2)):
        if p_res2[h] <= 0:
            p_res2[h] = 0

    return [p_res2, param_res]


[abe, banan] = fitting_peaks_double_interval(freq, p, 5)
print(banan)
print(banan[0])