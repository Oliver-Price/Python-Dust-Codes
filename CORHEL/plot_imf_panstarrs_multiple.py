# -*- coding: utf-8 -*-
# 1 r         2 lon       3 lat       4 N         5 T         6 V_r       7 V_lat     8 V_lon     9 B_r       10 B_lat    11 B_lon    12 polB     13 N*r^2    14 P        15 P*r^2     
# AU          deg         deg         cm^-3       K           km/s        km/s        km/s        nT          nT          nT          []          AU^2cm^-3   nPa         ???

import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

dir0 = r"C:\PhD\Python\Python-Dust-Codes\CORHEL\basedata"
dir1 = r"C:\PhD\Python\Python-Dust-Codes\CORHEL\dust_1.7_data"
dir2 = r"C:\PhD\Python\Python-Dust-Codes\CORHEL\dust_0.8_data"

data0 = np.empty((0,16))
data1 = np.empty((0,16))
data2 = np.empty((0,16))

for j in os.listdir(dir0):
    k = os.path.join(dir0,j)
    csv_temp = np.genfromtxt (k, delimiter=",")
    csv_temp = csv_temp[~np.isnan(csv_temp).all(axis=1)]
    data0 = np.append(data0,csv_temp, axis=0)
    
for j in os.listdir(dir1):
    k = os.path.join(dir1,j)
    csv_temp = np.genfromtxt (k, delimiter=",")
    csv_temp = csv_temp[~np.isnan(csv_temp).all(axis=1)]
    data1 = np.append(data1,csv_temp, axis=0)
    
for j in os.listdir(dir2):
    k = os.path.join(dir2,j)
    csv_temp = np.genfromtxt (k, delimiter=",")
    csv_temp = csv_temp[~np.isnan(csv_temp).all(axis=1)]
    data2 = np.append(data2,csv_temp, axis=0)   

times0 = Time(data0[:,0],format='jd')
taxis0 = (times0.plot_date)

times1 = Time(data1[:,0],format='jd')
taxis1 = (times1.plot_date)

times2 = Time(data2[:,0],format='jd')
taxis2 = (times2.plot_date)

t_start = Time('2013-03-11T00:00',format='isot').plot_date
t_end = Time('2013-03-17T00:00',format='isot').plot_date


#%%
#Br Bphi Btheta Vr Vphi Vtheta N T
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator) 

f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, sharex=True, figsize=(7,10))

#ax1.set_title("tstart + " to " + tend)
ax1.set_xlim(t_start, t_end)
f.autofmt_xdate()  

ax1.axhline(y=0,color='lightgrey',linestyle='dotted')
ax1.plot_date(taxis1,data1[:,9],'magenta',linestyle='--',lw=3)
ax1.plot_date(taxis2,data2[:,9],'lime',linestyle='--',lw=3)
ax1.plot_date(taxis0,data0[:,9],'k-',lw=3)

ax1.set_ylim(-25,25)
ax1.set_ylabel(r"$B_{r}$ (nT)") 

ax2.axhline(y=0,color='lightgrey',linestyle='dotted')
ax2.plot_date(taxis1,data1[:,10],'magenta',linestyle='--',lw=3)
ax2.plot_date(taxis2,data2[:,10],'lime',linestyle='--',lw=3)
ax2.plot_date(taxis0,data0[:,10],'k-',lw=3)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_ylim(-0.5,0.5)
ax2.set_ylabel(r"$B_{\phi}$ (nT)") 

ax3.axhline(y=0,color='lightgrey',linestyle='dotted')
ax3.plot_date(taxis1,data1[:,11],'magenta',linestyle='--',lw=3)
ax3.plot_date(taxis2,data2[:,11],'lime',linestyle='--',lw=3)
ax3.plot_date(taxis0,data0[:,11],'k-',lw=3)
ax3.set_ylim(-10,10)
ax3.set_ylabel(r"$B_{\theta}$ (nT)") 

ax4.plot_date(taxis1,data1[:,6],'magenta',linestyle='--',lw=3)
ax4.plot_date(taxis2,data2[:,6],'lime',linestyle='--',lw=3)
ax4.plot_date(taxis0,data0[:,6],'k-',lw=3)
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.set_ylim(100,700)
ax4.yaxis.set_major_locator(MultipleLocator(200))
ax4.set_ylabel(r"$V_{r}$ (km/s)") 

ax5.axhline(y=0,color='lightgrey',linestyle='dotted')
ax5.plot_date(taxis1,data1[:,7],'magenta',linestyle='--',lw=3)
ax5.plot_date(taxis2,data2[:,7],'lime',linestyle='--',lw=3)
ax5.plot_date(taxis0,data0[:,7],'k-',lw=3)
ax5.set_ylim(-10,10)
ax5.set_ylabel(r"$V_{\phi}$ (km/s)") 

ax6.axhline(y=0,color='lightgrey',linestyle='dotted')
ax6.plot_date(taxis1,data1[:,8],'magenta',linestyle='--',lw=3)
ax6.plot_date(taxis2,data2[:,8],'lime',linestyle='--',lw=3)
ax6.plot_date(taxis0,data0[:,8],'k-',lw=3)
ax6.yaxis.tick_right()
ax6.yaxis.set_label_position("right")
ax6.set_ylim(-30,30)
ax6.set_ylabel(r"$V_{\theta}$ (km/s)") 

ax7.plot_date(taxis1,data1[:,4],'magenta',linestyle='--',lw=3)
ax7.plot_date(taxis2,data2[:,4],'lime',linestyle='--',lw=3)
ax7.plot_date(taxis0,data0[:,4],'k-',lw=3)
ax7.set_ylim(-20,220)
ax7.yaxis.set_major_locator(MultipleLocator(100))
ax7.set_ylabel(r"$N$ ($cm^{-3}$)")

ax8.plot_date(taxis1,data1[:,5],'magenta',linestyle='--',lw=3)
ax8.plot_date(taxis2,data2[:,5],'lime',linestyle='--',lw=3)
ax8.plot_date(taxis0,data0[:,5],'k-',lw=3)
ax8.yaxis.tick_right()
ax8.yaxis.set_label_position("right")
ax8.set_ylim(-20,200000)
ax8.yaxis.set_major_locator(MultipleLocator(100000))
ax8.set_ylabel(r"T (K)") 

f.subplots_adjust(hspace=0.1)
f.subplots_adjust(left=0.11,right=0.88,top=0.95,bottom=0.08)
ax1.xaxis.set_minor_locator(MultipleLocator(1/4))
[a.tick_params(which='minor', color='gray') for a in f.axes]

plt.savefig('C:\PhD\Python\Python-Dust-Codes\CORHEL\Panstarrs_SW_All.png')