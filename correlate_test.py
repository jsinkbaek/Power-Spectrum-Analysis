"""
Script used to test the function as_f.autocorr, by passing it a generated spectrum, that has a known autocorrelation.
"""

import as_f
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 26})

# # Define test spectrum
x = np.linspace(0, 3000, 90000)
y = [np.sinc(0.1*(i-500))**2 + np.sinc(0.1*(i-1000))**2 + np.sinc(0.1*(i-1300))**2 for i in x]

# # Plot test spectrum
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend(['f(x) = sinc(0.1*(x-500))^2 + sinc(0.1*(x-1000))^2 + sinc(0.1*(x-1300))^2'], fontsize=12)
plt.show()

# # Perform autocorrelation
correlation = as_f.autocorr(y)

# # Plot autocorrelation
plt.plot(x, correlation)
plt.xlabel('displacement of x')
plt.ylabel('Autocorrelation')
plt.legend(['Autocorrellation of f as a function of displacement of x'])
plt.show()


