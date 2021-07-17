# -*- coding: utf-8 -*-
# 1 r         2 lon       3 lat       4 N         5 T         6 V_r       7 V_lat     8 V_lon     9 B_r       10 B_lat    11 B_lon    12 polB     13 N*r^2    14 P        15 P*r^2     
# AU          deg         deg         cm^-3       K           km/s        km/s        km/s        nT          nT          nT          []          AU^2cm^-3   nPa         ???

import os
#os.chdir(r"C:\PhD\Python\Python-Dust-Codes\CORHEL\basedata")
os.chdir(r"C:\PhD\Python\Python-Dust-Codes\CORHEL\dust_0.8_data")
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

data = np.empty((0,16))

for j in os.listdir():
    csv_temp = np.genfromtxt (j, delimiter=",")
    csv_temp = csv_temp[~np.isnan(csv_temp).all(axis=1)]
    data = np.append(data,csv_temp, axis=0)

data = data[2:-3,:]
times = Time(data[:,0],format='jd')
taxis = (times.plot_date) #-719163

#%%
'''
plt.axhline(y=0,color='lightgrey',linestyle='dotted')
plt.plot_date(taxis,data[:,9],'k-')
plt.xlim([taxis[0], taxis[-1]])

plt.axhline(y=0,color='lightgrey',linestyle='dotted')
plt.plot(data[:,0],data[:,9],'k-')'
'''

#%%
#Br Bphi Btheta Vr Vphi Vtheta N T
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator) 

f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, sharex=True, figsize=(7,10))
tstart = times.isot[0][0:16].replace("T"," ")
tend = times.isot[-1][0:16].replace("T"," ")

#ax1.set_title("tstart + " to " + tend)
ax1.set_xlim(taxis[0], taxis[-1])
f.autofmt_xdate()  

ax1.axhline(y=0,color='lightgrey',linestyle='dotted')
ax1.plot_date(taxis,data[:,9],'k-')
ax1.set_ylim(-25,25)
ax1.set_ylabel(r"$B_{r}$ (nT)") 

ax2.axhline(y=0,color='lightgrey',linestyle='dotted')
ax2.plot_date(taxis,data[:,10],'k-')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_ylim(-0.5,0.5)
ax2.set_ylabel(r"$B_{\phi}$ (nT)") 

ax3.axhline(y=0,color='lightgrey',linestyle='dotted')
ax3.plot_date(taxis,data[:,11],'k-')
ax3.set_ylim(-10,10)
ax3.set_ylabel(r"$B_{\theta}$ (nT)") 

ax4.plot_date(taxis,data[:,6],'k-')
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.set_ylim(100,700)
ax4.yaxis.set_major_locator(MultipleLocator(200))
ax4.set_ylabel(r"$V_{r}$ (km/s)") 

ax5.axhline(y=0,color='lightgrey',linestyle='dotted')
ax5.plot_date(taxis,data[:,7],'k-')
ax5.set_ylim(-10,10)
ax5.set_ylabel(r"$V_{\phi}$ (km/s)") 

ax6.axhline(y=0,color='lightgrey',linestyle='dotted')
ax6.plot_date(taxis,data[:,8],'k-')
ax6.yaxis.tick_right()
ax6.yaxis.set_label_position("right")
ax6.set_ylim(-30,30)
ax6.set_ylabel(r"$V_{\theta}$ (km/s)") 

ax7.plot_date(taxis,data[:,4],'k-')
ax7.set_ylim(-20,220)
ax7.yaxis.set_major_locator(MultipleLocator(100))
ax7.set_ylabel(r"$N$ ($cm^{-3}$)")

ax8.plot_date(taxis,data[:,5],'k-')
ax8.yaxis.tick_right()
ax8.yaxis.set_label_position("right")
ax8.set_ylim(-20,200000)
ax8.yaxis.set_major_locator(MultipleLocator(100000))
ax8.set_ylabel(r"T (K)") 

f.subplots_adjust(hspace=0.1)
f.subplots_adjust(left=0.15,right=0.85,top=0.95,bottom=0.08)
ax1.xaxis.set_minor_locator(MultipleLocator(1/4))
[a.tick_params(which='minor', color='gray') for a in f.axes]

plt.savefig('C:\PhD\Python\Python-Dust-Codes\CORHEL\PanstarrsSWdust0p8.png')