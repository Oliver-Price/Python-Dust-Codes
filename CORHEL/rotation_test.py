# -*- coding: utf-8 -*-
import os
os.chdir("C:\PhD\Python\Python-Dust-Codes\CORHEL")
import numpy as np
import matplotlib.pyplot as plt

csv9 = np.genfromtxt ('bpol-2013-3-12-09-40-18.csv', delimiter=",")
csv15 = np.genfromtxt ('bpol-2013-3-12-15-40-18.csv', delimiter=",")

lon9 = np.resize(csv9[:,1],(59,180))
lat9 = np.resize(csv9[:,2],(59,180))
pol9 = np.resize(csv9[:,3],(59,180))

lon15 = np.resize(csv15[:,1],(59,180))
lat15 = np.resize(csv15[:,2],(59,180))
pol15 = np.resize(csv15[:,3],(59,180))

cor = 360/27.2753/24

fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)
axs[0].pcolormesh(lon9+cor*6,lat9,pol9)
axs[1].pcolormesh(lon15,lat15,pol15)

