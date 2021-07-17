# -*- coding: utf-8 -*-
from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection="3d")

ax.plot3D(comveceq[:,6],comveceq[:,7],comveceq[:,8],c='b')
ax.plot3D(obsveceq[:,6],obsveceq[:,7],obsveceq[:,8],c='g')
ax.scatter3D(obsveceq[-1,6],obsveceq[-1,7],obsveceq[-1,8],c='g')
ax.scatter3D(comveceq[-1,6],comveceq[-1,7],comveceq[-1,8],c='b')
ax.scatter3D(0,0,0, c='orange')

for a in range(np.shape(simres)[1]):
    locs = np.where(simres[:,a,14])
    ax.plot3D(simres[locs,a,4].ravel(),simres[locs,a,5].ravel(),simres[locs,a,6].ravel(),c='r')
    
plt.show()

pos = np.empty_like(simres[:,:,4:7])
pos[...,0] = simres[...,4] - obsveceq[int(simres[0,0,3]),6]
pos[...,1] = simres[...,5] - obsveceq[int(simres[0,0,3]),7]
pos[...,2] = simres[...,6] - obsveceq[int(simres[0,0,3]),8]

rc = np.linalg.norm(pos,axis=2)#*simres[...,14]
scale = rc*(70/60/60*np.pi/180*1.5e8)