# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

bvalsav = r'C:\PhD\Python\fake_directory\bvals.npy'
tvalsav = r'C:\PhD\Python\fake_directory\tvals.npy'
srcolsav = r'C:\PhD\Python\fake_directory\srcolors.npy'

srcolors = np.load(srcolsav)
bvals = np.load(bvalsav)
tvals = np.load(tvalsav)

sr_long = np.empty((1000000,3), dtype=float)

a = 0
for t in range(0,1000):
    for b in range(0,1000):
        sr_long[a,0] = tvals[t]
        sr_long[a,1] = bvals[b]
        sr_long[a,2] = srcolors[t,b] 
        a+=1
        
#%%
        
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
tes = ax.scatter(sr_long[:,0],sr_long[:,1],c=sr_long[:,2],s=20, lw = 0)
        