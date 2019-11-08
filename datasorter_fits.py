from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
from ps_f import writer2col

# # Indicate and show data structure
fits_file = 'Datafiles/GOLF_22y_MEAN.fits'
hdu_list = fits.open(fits_file)
hdr = hdu_list


# # Read file and save as arrays
with fits.open(fits_file, mode="readonly") as hdulist:
    tess_bjd = hdulist[1].data['TIME']
    sap_flux = hdulist[1].data['SAP_FLUX']
    pdcsap_flux = hdulist[1].data['PDCSAP_FLUX']
    qual_flags = hdulist[1].data['QUALITY']

# Locate quality flags greater than 2
where_qgt0 = np.where(qual_flags > 2)[0]
where_not_qgt0 = np.where(qual_flags <= 2)[0]


# # Plot

# Start figure and axis
fig, ax = plt.subplots()

# Plot corrected (PDC) timeseries
ax.plot(tess_bjd, pdcsap_flux, 'k.')

# Overplot fluxes with quality flags greater than 2
ax.plot(tess_bjd[where_qgt0], pdcsap_flux[where_qgt0], 'ro')

# Label axis
ax.set_ylabel('PDCSAP Flux (e-/s)')
ax.set_xlabel('Time (TBJD)')

# Show
plt.show()


# # Limit data to fluxes with quality flag less or equal to 2
[tess_bjd, sap_flux, pdcsap_flux] = [tess_bjd[where_not_qgt0], sap_flux[where_not_qgt0],
                                     pdcsap_flux[where_not_qgt0]]


# # Find aN (not NaN) values and remove NaN
where_an = np.where(~np.isnan(tess_bjd))[0]
[tess_bjd, sap_flux, pdcsap_flux] = [tess_bjd[where_an], sap_flux[where_an], pdcsap_flux[where_an]]


# # Find median value, divide and center
pdcsap_median = np.median(pdcsap_flux)
pdcsap_flux[:] = [x / pdcsap_median - 1 for x in pdcsap_flux]


# # Plot to check
plt.plot(tess_bjd, pdcsap_flux, 'b.')
plt.show()

# # Save as .dat file
inpt = input('What name for the saved data? (exclude .dat)')
writer2col(inpt, tess_bjd, pdcsap_flux)
