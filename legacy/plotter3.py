import matplotlib.pyplot as plt
import ps_f
import numpy as np
plt.rcParams.update({'font.size': 30})

# # Load data
data1 = ps_f.preader('filtered_star2_extra_3std.dat')
data2 = ps_f.preader('filtered_star2_extra_old.dat')

# # Save in lists

# Filtered
t1 = []
c1 = []  # corrected flux ppm
m1 = []  # median corrected flux
v1 = []  # corrected flux variance from median
for row in data1:
    t1.append(row[0])
    c1.append(row[1])
    m1.append(row[2])
    v1.append(row[3])

# Not filtered (old)
t2 = []
c2 = []
m2 = []
v2 = []
for row in data2:
    t2.append(row[0])
    c2.append(row[1])
    m2.append(row[2])
    v2.append(row[3])

# Convert to days
t1[:] = [x / 86400 for x in t1]
t2[:] = [x / 86400 for x in t2]

f = plt.figure()

ax01 = f.add_subplot(1, 2, 1)
ax01.set_ylabel('Variance')
plt.xticks([])
plt.yticks([])
ax02 = f.add_subplot(1, 2, 2)
ax02.yaxis.set_label_position("right")
plt.xticks([])
plt.yticks([])
ax02.set_ylabel('Sorted Variance')

ax1 = f.add_subplot(2, 2, 1)
plt.plot(t2, v2, 'r.')
ax1.yaxis.tick_right()
ax1.set_yticklabels([])

ax2 = f.add_subplot(2, 2, 2)
plt.plot(sorted(v2), 'r.')

ax3 = f.add_subplot(2, 2, 3, sharey=ax1)
plt.xlabel('Time (d)')
plt.plot(t1, v1, 'r.')
ax3.yaxis.tick_right()
ax3.set_yticklabels([])

ax4 = f.add_subplot(2, 2, 4, sharey=ax2)
plt.xlabel('Element')

plt.plot(sorted(v1), 'r.')



plt.show()
