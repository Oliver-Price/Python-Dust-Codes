# -*- coding: utf-8 -*-
#program to plot orbits of various spacecraft
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib

keys = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
values = range(1,13)
mondict = dict(zip(keys, values))

def loadorbit(datapath):

    with open(datapath, "r") as dat:
        
        #reads raw data as string 
        rawdata = dat.readlines()
        datasize = int(len(rawdata)/2)
        data = np.empty((datasize,9),dtype = float)
        
        #converts string to numpy array
        for hrow in range(0,datasize,1):
            data[hrow,0] = rawdata[2*hrow][0:17]            #jd
            data[hrow,1] = rawdata[2*hrow][25:29]           #year
            data[hrow,2] = mondict[rawdata[2*hrow][30:33]]  #month
            data[hrow,3] = rawdata[2*hrow][34:36]           #day
            data[hrow,4] = rawdata[2*hrow][37:39]           #hour
            data[hrow,5] = rawdata[2*hrow][40:42]           #minute
            data[hrow,6] = rawdata[(hrow*2)+1][4:26]
            data[hrow,7] = rawdata[(hrow*2)+1][30:52]
            data[hrow,8] = rawdata[(hrow*2)+1][56:78]
        
    return data

filelocs = 'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Orbit_Plots'
earthloc = os.path.join(filelocs,'EARTH_ORBIT.txt')
mcnaughtloc = os.path.join(filelocs,'MCNAUGHT_ORBIT.txt')
soholoc = os.path.join(filelocs,'SOHO_ORBIT.txt')
steraloc = os.path.join(filelocs,'STEREO_A_ORBIT.txt')
sterbloc = os.path.join(filelocs,'STEREO_B_ORBIT.txt')

earth_data = loadorbit(earthloc)
mcnaught_data = loadorbit(mcnaughtloc)
soho_data = loadorbit(soholoc)
stereo_a_data = loadorbit(steraloc)
stereo_b_data = loadorbit(sterbloc)

plt.figure(figsize=(5,5))
tes1 = plt.plot(earth_data[:,6], earth_data[:,7])
tes2 = plt.plot(mcnaught_data[:,6], mcnaught_data[:,7])
tes3 = plt.plot(soho_data[:,6], soho_data[:,7])
tes4 = plt.plot(stereo_a_data[:,6], stereo_a_data[:,7])
tes5 = plt.plot(stereo_b_data[:,6], stereo_b_data[:,7])
axes = plt.gca()
axes.set_xlim([-1,1])
axes.set_ylim([-1,1])
