import ps_f
import matplotlib.pyplot as plt
import numpy as np
import as_f
from pandas import Series

plt.rcParams.update({'font.size': 30})  # Changes font size for all plots

# # Data load from file
# Input data file name
print('Input data to load, exclude extension (.dat)')
#inpt1 = input()
inpt1 = 'md_sn'
# Call reader function to pull data into list
data = ps_f.dreader(inpt1 + '.dat')

# Preload data lists
freq = []
p = []

# Append data to lists
for row in data:
    freq.append(row[0])
    p.append(row[1])
muHz = 0.000001
freq_muHz = [x / muHz for x in freq]

as_f.smoother(freq_muHz, p)
