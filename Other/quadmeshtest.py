# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

nx, ny = (100,100)
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
xv, yv = np.meshgrid(x, y)

c = xv-yv

plt.pcolormesh(xv,yv,c)
plt.xlabel("top")
plt.ylabel("side")