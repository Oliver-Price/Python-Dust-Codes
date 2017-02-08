import matplotlib.pyplot as plt
import numpy as np

Fs = 10
sample = 1000
x = np.arange(sample)
y = x/(np.sqrt(10+100*x**2))

#y = np.sin(2 * np.pi * f * x / Fs)
plt.plot(x, y)
plt.xlabel('voltage(V)')
plt.ylabel('sample(n)')
plt.yscale('log')
plt.xscale('log')
plt.show()