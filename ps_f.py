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


def writer(fname, *args):
    fname = os.path.join('Datafiles', fname + '.dat')
    np.savetxt(fname, np.c_[args])


def reader(fname):
    fname = os.path.join('Datafiles', fname + '.dat')
    arr = np.loadtxt(fname)
    return np.array_split(arr, len(arr[0, :]), axis=1)


