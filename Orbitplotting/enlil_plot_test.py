# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime
import os

folder = r"C:\PhD\Python\Python-Dust-Codes\Orbitplotting"
sys.path.append(folder)

from coordinate_transform_test import *

filein = os.path.join(folder,'pointdata_137825026838_20130313_033852.csv')

data = np.genfromtxt(filein,delimiter=',')
y = -data[:,0]*np.sin(np.radians(data[:,1]))
x = -data[:,0]*np.cos(np.radians(data[:,1]))
#plt.scatter(x,y,c=data[:,3],cmap='jet')

cyear = 2013
cmonth = 3
cday = 13
chour = 3
cmin = 38
csec = 52

ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                            chour , cmin, csec))  
    
[x_heeq,y_heeq,z_heeq] = heeq2heeqxyz(data[:,0],data[:,2],data[:,1])

n = get_n(ctime)
L = get_L(n)
g = get_g(n)
lmb = get_lambda(L,g)
omg = get_omega(ctime)
tht = get_theta(lmb,omg)

[x_hae,y_hae,z_hae] = heeq2hae(x_heeq,y_heeq,z_heeq,tht,omg)

#%%

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_hae,y_hae,zs=z_hae,c=data[:,3])
