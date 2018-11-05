# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

nx, ny = (100,100)
x = np.linspace(-2, 2, nx)
y = np.linspace(-1, 1, ny)
xv, yv = np.meshgrid(x, y)

c = xv-yv

plt.pcolormesh(xv,yv,c)
plt.xlabel("top")
plt.ylabel("side")

#%%

w = 5
h = 5

im_np = np.random.rand(h, w)

fig = plt.figure(frameon=False)
fig.set_size_inches(w,h)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(im_np)
fig.savefig(r'C:\PhD\Python\Python-Dust-Codes\Other\figure.png', dpi=100)