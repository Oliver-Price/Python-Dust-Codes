# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime
import os
import astropy.time
import time

folder = r"C:\PhD\Python\Python-Dust-Codes\Orbitplotting"
sys.path.append(folder)

from coordinate_transform_test import *

#earth
filein = os.path.join(folder,'pointdata_211205033915.csv')

#stereo-a
#filein = os.path.join(folder,'pointdata_732405045063.csv')

earthdata = np.genfromtxt(filein,delimiter=',')

cyear = 2013
cmonth = 2
cday = 22
chour = 3
cmin = 0
csec = 0

ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                            chour , cmin, csec))

stime_arr = np.full(earthdata.shape[0],ctime.jd)
stime = astropy.time.Time(stime_arr,format='jd')
dtime = astropy.time.TimeDelta(earthdata[:,0],format='jd')

times = stime + dtime

[x_heeq,y_heeq,z_heeq] = rlatlon2xyz(earthdata[:,1],earthdata[:,3],earthdata[:,2])

n = get_n(times)
L = get_L(n)
g = get_g(n)
lmb = get_lambda(L,g)
omg = get_omega(times)
tht = get_theta(lmb,omg)

[x_hae,y_hae,z_hae] = heeq2hae(x_heeq,y_heeq,z_heeq,tht,omg)

#%%
z = 300#tht - 90
    
cosz = np.cos(np.radians(z))
sinz = np.sin(np.radians(z))

x_corr = x_hae*cosz - y_hae*sinz
y_corr = x_hae*sinz + y_hae*cosz

plt.plot(times.jd,x_corr,'r')
plt.plot(times.jd,y_corr,'g')
plt.plot(times.jd,z_hae,'b')

#%%
plt.plot(obsveceq[:,0],obsveceq[:,6],'r--')
plt.plot(obsveceq[:,0],obsveceq[:,7],'g--')
plt.plot(obsveceq[:,0],obsveceq[:,8],'b--')

#%%
'''
[r,lat,lon] = xyz2rlatlon(x_hae,y_hae,z_hae)

plt.plot(times.jd,r,'r')
plt.plot(times.jd,lat,'g')
plt.plot(times.jd,lon,'b')
'''