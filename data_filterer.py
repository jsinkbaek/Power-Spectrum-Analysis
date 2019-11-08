"""
Script to be used for filtering of corrected flux data files generated by datasorter.py.
- Loads corrected flux and time from file.
- For each data-point, calculates a median flux value from the surrounding 9 points (including itself, 4 to each side).
- Plots variance of data from median value.
- Optional: Allows you to graphically select a specific region to exclude from the data-set.
- Optional: Can automatically delete data that is surrounded by bad data points. Not recommended, not fully tested.
- Calculates standard deviation of variance from median.
- After this, will automatically delete data points that has a variance larger than n standard deviations (n chosen by
  text prompt).
- Saves filtered data to file.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 22})  # Changes font size on all plots


def reader(filename):  # Assumes only two columns, time and cflux, function that reads data from file.
    A = []
    f = open(os.path.join('Datafiles', filename), 'r')
    for line in f:
        columns = line.split()
        columns[0] = float(columns[0])
        columns[1] = float(columns[1])
        A.append(columns)
    f.close()
    return A


def writer(filename, time, cflux, mcflux, variance):  # Function that writes the results to an indicated file
    open(os.path.join('Datafiles', filename + '.dat'), 'w').close()
    k = len(time)
    for i in range(0, k):
        content = str(time[i]) + '   ' + str(cflux[i]) + '   ' + str(mcflux[i]) + '   ' + str(variance[i]) + '\n'
        with open(os.path.join('Datafiles', filename + '.dat'), 'a') as f:
            f.write(content)


print('What is the filename (without extension)?')
inpt1 = input()  # Input for data reader
B = reader(inpt1 + '.dat')  # reader function call to load data


time = []
cflux = []
for row in B:  # Loop that saves data in the variables tid and cflux
    time.append(row[0])
    cflux.append(row[1])
start = time[0]
time[:] = [x - start for x in time]  # starting point correction (if different than 0)
time_days = [x for x in time]
time[:] = [x * 86400 for x in time]  # change to time in seconds


# # Sort out data points at ends, not necessary but there is a higher chance of bad data here from the satellite
u = 25
del time[-u:-1]
del cflux[-u:-1]
del time[0:u]
del cflux[0:u]


# # Median points
n = 5  # Decides how many points are grouped together. Can be changed if necessary
mcflux = [0] * len(cflux)  # preload variable
b = True  # check used to see if the loop has reached the end
x = n  # Start point
while b is True:  # Loop to create median points from surrounding 8 points
    mcflux[x] = np.median(cflux[x-n+1:x+n])
    if x == len(cflux)-n:
        b = False
    else:
        x += 1
del mcflux[-n:-1]  # Not pretty, but done to make sure that right data is used
del mcflux[0:n]
tid_mc = time[n:-n]  # New shorter variable matching in length with mcflux
del mcflux[-1]


# # Median plot
plt.plot(time, cflux, 'b.')
plt.plot(time, cflux, 'b--', linewidth=0.5)
plt.xlabel('Days from t0')
plt.ylabel('Corrected flux ppm')
plt.plot(tid_mc, mcflux, 'g.')
plt.legend(['Data', 'Median of data from surrounding points'])

plt.show()


# #  Distribution of points
cflux_short = cflux  # New shorter variable matching in length with mcflux
del cflux_short[-n:-1]
del cflux_short[0:n]
del cflux_short[-1]

variance = [i - j for i, j in zip(cflux_short, mcflux)]  # subtracting mcflux from cflux_short to find median variance
h = sorted(variance)  # Variance sorted from low to high value


# # Plot variance
m = len(variance)
plt.plot(tid_mc, variance, 'r.')
plt.xlabel('time (s)')
plt.ylabel('Variance')
plt.show()

plt.plot(sorted(variance), 'r.')
plt.ylabel('sorted variance')
plt.show()


# # Deletion of data points by manual selection of interval
print('Do you want to remove points in a region manually? (y/n)')
inpt4 = input()
if inpt4 == 'y':
    # # Graphical input selection of region to delete
    plt.plot(variance, 'b.')
    plt.xlabel('List element')
    plt.ylabel('Variance')
    plt.title('Select region which will be deleted (only 1 region, two points)')
    coords = plt.ginput(n=2, timeout=0, show_clicks=True, mouse_add=1, mouse_stop=3, mouse_pop=2)  # graphical input

    # # Conversion of from input coordinates to x-axis integer coordinates, x-axis is list index
    [coord1, coord2] = [coords[0], coords[1]]
    coord1 = int(coord1[0])
    coord2 = int(coord2[0])

    del_select = list(range(coord1, coord2))  # List with deletion indices to be used for plotting

    # # Plot selected region
    plt.plot(variance, 'r.')
    plt.plot(del_select, variance[coord1:coord2], 'g.')
    plt.show()

    # # Delete selected region from data set
    del variance[coord1:coord2]
    del tid_mc[coord1:coord2]
    del cflux_short[coord1:coord2]
    del mcflux[coord1:coord2]


# # Deletion of data points by mark-filtering (Experimental, not thoroughly tested)
if inpt4 != 'y':
    print('Do you want to remove points by automatic marking (y/n)? Warning: Experimental, not tested thoroughly.')
    inpt5 = input()

    if inpt5 == 'y':
        std_limit = 3*np.std(h)  # Sets the lowest variance necessary to result in a mark
        marks_limit = 3  # Sets the amount of marks necessary for the script to delete the interval
        k = len(variance) - 2*n - 2  # Current step, n is about half the amount of values to include in the interval
        while k >= 0:  # Loop that goes through the whole data set (from end to start)
            delete = False
            marks = 0
            temp = variance[k-2*n:k+2*n]  # the interval to be checked
            for j in range(0, len(temp)):  # Checks and marks the interval depending on std_limit
                if temp[j] <= -std_limit or temp[j] >= std_limit:
                    marks += 1
                if marks >= marks_limit:  # Checks if the interval has more marks than acceptable
                    delete = True
            if delete is True:  # Deletes elements in the interval from data sets if marked for deletion
                del variance[k-2*n:k+2*n]
                del tid_mc[k-2*n:k+2*n]
                del cflux_short[k-2*n:k+2*n]
                del mcflux[k-2*n:k+2*n]
                k = k - 4 * n  # Changes step to accurately reflect deletion
            elif delete is False:
                k = k - 1  # Changes step by one if no deletion was done
            else:
                print('Error')  # Unknown error indicator

# # Plotting of variance after possible deletion of data
plt.plot(tid_mc, variance, 'r.')
plt.xlabel('time (s)')
plt.ylabel('Variance')
plt.show()

# # Plotting of sorted variance after possible deletion of data
plt.plot(sorted(variance), 'r.')
plt.ylabel('sorted variance')
plt.show()


# # Input to decide what the upper limit for deleting a point is (indicated in stand dev, usually 3 or 4 is used)
print('How many standard deviations away do we short out data?')
inpt2 = float(input())  # Must be an integer


# #  Deletion of points further than inpt2 std away in variance
stdev = np.std(variance)  # Standard deviation of variance
m = len(variance)

variance_old = [x for x in variance]
tid_mc_old = [x for x in tid_mc]
cflux_short_old = [x for x in cflux_short]
mcflux_old = [x for x in mcflux]

for i in range(1, m+1):  # Checks and deletes if each point has a higher numerical value than the indicated limit
    if abs(variance[m-i]) > inpt2*stdev:
        del variance[m-i]
        del cflux_short[m-i]
        del mcflux[m-i]
        del tid_mc[m-i]

# # Plot new sorted variance
plt.plot(sorted(variance), 'r.')
plt.ylabel('sorted variance after discarding points')
plt.show()

# # Plot new variance vs. time
plt.plot(tid_mc, variance, 'r.')
plt.xlabel('time (s)')
plt.ylabel('Variance')
plt.show()

# # Plot new data set (cflux)
plt.plot(tid_mc, cflux_short, 'b.')
plt.plot(tid_mc, mcflux, 'g.')
plt.xlabel('Time (s)')
plt.ylabel('Corrected flux (ppm)')
plt.show()

# # Writer call to save data
print('What is the filename (without extension)?')
inpt3 = input()
writer(inpt3, tid_mc, cflux_short, mcflux, variance)
inpt4 = input('What is the filename of old data')
writer(inpt4, tid_mc_old, cflux_short_old, mcflux_old, variance_old)
