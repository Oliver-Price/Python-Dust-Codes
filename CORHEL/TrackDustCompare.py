# -*- coding: utf-8 -*-
from astropy.time import Time
import easygui
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\Orbitplotting")

from orbitdata_loading_functions import *
from FP_plot_functions import *
from FP_diagnostics import *
from imagetime_methods import *
from conversion_routines import *
from io_methods import *
from simulation_setup import *
from particle_sim import *
from enlil_cooordinate_transform import xyz2rlatlon
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import astropy.units as u

t_dust1 = Time('2013-03-05T12:00',format='isot')
beta1 = 0.8

t_dust2 = Time('2013-03-08T04:53',format='isot')
beta2 = 1.7


#choosing comet data to use
inputfile = r'C:\\PhD\\Comet_data\\Input_files\\Input file_PANSTARRS_c2011l4_pt1.txt'

#reading main comet parameters
with open(inputfile, "r") as c:
    cdata = c.readlines()
    comname = cdata[30][12:]
    comdenom = cdata[31][13:-2]
    imagedir = cdata[24][18:-2]
    orbitdir = cdata[25][23:-2]
    idlsav = cdata[26][25:-2]
    pysav = cdata[27][24:-2]
    obslocstr = cdata[34][19:]
    horiztag = cdata[40][10:]

obsloc = 'Stereo_B'
 
#import the orbit data
comveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'eq')

comcell1 = abs(comveceq[:,0] - t_dust1.jd).argmin()
comjds = Time(comveceq[:,0],format='jd').squeeze()

simt = 20*60*24 #10 days in minutes
nperday = 48
tstart1 = comveceq[comcell1,6:12]

(times1,positions1) = part_sim_fine_track(beta1, simt, nperday, tstart1)
jds1 = Time(t_dust1.jd+times1/1440,format='jd').squeeze()

comcell2 = abs(comveceq[:,0] - t_dust2.jd).argmin()

simt = 20*60*24 #10 days in minutes
nperday = 48
tstart2 = comveceq[comcell2,6:12]

(times2,positions2) = part_sim_fine_track(beta2, simt, nperday, tstart2)
jds2 = Time(t_dust2.jd+times2/1440,format='jd').squeeze()

#%%
fig = plt.figure()
ax = plt.axes(projection="3d")

ax.plot3D(comveceq[:,6], comveceq[:,7], comveceq[:,8], 'black')
ax.plot3D(positions1[:,0], positions1[:,1], positions1[:,2], 'green')
ax.plot3D(positions2[:,0], positions2[:,1], positions2[:,2], 'magenta')