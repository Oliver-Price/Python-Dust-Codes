# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy

# Number of samplepoints
N = 800
# sample spacing
T = 1.0 / N
x = np.linspace(0.0, N*T, N)
y = np.sin(50*2*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig, ax = plt.subplots()
ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
plt.show()
