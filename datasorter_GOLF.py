from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
from ps_f import writer2col

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
print('time_jd.shape', time_jd.shape)
print('flux_data.shape', flux_data.shape)


inpt = input('What name for the saved data? (exclude .dat)')
print(int(995*len(flux_data)/1000))
writer2col(inpt, time_jd[int(995*len(time_jd)/1000):len(time_jd)], flux_data[int(995*len(flux_data)/1000):len(flux_data)])
