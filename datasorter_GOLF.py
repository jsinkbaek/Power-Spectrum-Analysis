from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
from ps_f import writer

# # Indicate and show data structure
fits_file = 'Datafiles/GOLF_22y_MEAN.fits'
hdu_list = fits.open(fits_file)
hdr = hdu_list[0].header
hdu_list.info()
print(repr(hdr))


# # Read file and save as arrays
jdstart = hdr['JDSTART']
jdend = hdr['JDEND']
cadence = hdr['CADENCE']
flux_data = hdu_list[0].data

time_jd = np.linspace(jdstart, jdend, (jdend-jdstart)/(cadence/86400))

inpt = input('What name for the saved data? (exclude .dat)')
k = 970/1000

writer(inpt, time_jd[int(k*len(time_jd)):len(time_jd)], flux_data[int(k*len(flux_data)):len(flux_data)])
