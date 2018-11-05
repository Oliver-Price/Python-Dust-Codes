#*************************************************************
#Program to visualise fits image of comet and overlay a
#finson-probstein diagram according to it's orbital parameters
#*************************************************************

import easygui
import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
import astropy.time
import pickle
import matplotlib.path as mplPath
from PIL import Image, ImageDraw, ImageFont
import webbrowser
import matplotlib.pyplot as plt

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from orbitdata_loading_functions import *
from FP_plot_functions import *
from FP_diagnostics import *
from imagetime_methods import *
from conversion_routines import *
from io_methods import *
from simulation_setup import *
from particle_sim import *

#%%**********************
#FIRST CELL - GET DATA IN
#************************

#choosing comet data to use
inputfilefolder = r"C:\PhD\Comet_data\Input_files\*pt1.txt"
inputfile = easygui.fileopenbox(default = inputfilefolder)

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
    an = float(cdata[33][16:])
    inc = float(cdata[32][13:])

#choose observer locations
obsloc = 'Earth'

#import the orbit data
comvec = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag)

comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')

#%% equatorial plane to ecliptic

i = -23.43929
cosi = np.cos(np.radians(i))
sini = np.sin(np.radians(i))

uneq = np.empty_like(comveceq[:,6:12])

uneq[:,0] = comveceq[:,6]
uneq[:,1] = comveceq[:,7]*cosi - comveceq[:,8]*sini
uneq[:,2] = comveceq[:,7]*sini + comveceq[:,8]*cosi
uneq[:,3] = comveceq[:,9]
uneq[:,4] = comveceq[:,10]*cosi - comveceq[:,11]*sini
uneq[:,5] = comveceq[:,10]*sini + comveceq[:,11]*cosi   


#%% ecliptic to orbital plane
        
#z = -87.41489440044852
#x = 77.83726262785106
z = 180 - an
x = inc

cosz = np.cos(np.radians(z))
sinz = np.sin(np.radians(z))
cosx = np.cos(np.radians(x))
sinx = np.sin(np.radians(x))

temp = np.empty_like(comvec[:,6:9])
orbp = np.empty_like(comvec[:,6:9])

temp[:,0] = comvec[:,6]*cosz - comvec[:,7]*sinz
temp[:,1] = comvec[:,6]*sinz + comvec[:,7]*cosz
temp[:,2] = comvec[:,8]

orbp[:,0] = temp[:,0]
orbp[:,1] = temp[:,1]*cosx - temp[:,2]*sinx
orbp[:,2] = temp[:,1]*sinx + temp[:,2]*cosx
 
        
#%%

cross_vec = np.cross(comvec[0,6:9],comvec[-1,6:9])
pn = cross_vec/np.linalg.norm(cross_vec)

#%% ecliptic to orbital plane speeds

tempv = np.empty_like(comvec[:,6:9])
orbv = np.empty_like(comvec[:,6:9])

tempv[:,0] = comvec[:,9]*cosz - comvec[:,10]*sinz
tempv[:,1] = comvec[:,9]*sinz + comvec[:,10]*cosz
tempv[:,2] = comvec[:,11]

orbv[:,0] = tempv[:,0]
orbv[:,1] = tempv[:,1]*cosx - tempv[:,2]*sinx
orbv[:,2] = tempv[:,1]*sinx + tempv[:,2]*cosx

#%%

theta_hat = np.empty_like(comvec[:,6:8])

rho = np.sqrt(orbp[:,0]**2 + orbp[:,1]**2)
theta_hat[:,0] = -orbp[:,1]/rho[:]
theta_hat[:,1] = orbp[:,0]/rho[:]

#%%

v_theta = orbv[:,0]*theta_hat[:,0] + orbv[:,1]*theta_hat[:,1]

z = 180 - an
x = inc

#%% equatorial to orbital plane

z = 180 - an
x = inc

cosz = np.cos(np.radians(z))
sinz = np.sin(np.radians(z))
cosx = np.cos(np.radians(x))
sinx = np.sin(np.radians(x))

temp1 = np.empty_like(comveceq[:,6:9])
temp2 = np.empty_like(comveceq[:,6:9])
orbp2 = np.empty_like(comveceq[:,6:9])

temp1[:,0] = comveceq[:,6]
temp1[:,1] = comveceq[:,7]*cosi - comveceq[:,8]*sini
temp1[:,2] = comveceq[:,7]*sini + comveceq[:,8]*cosi   

temp2[:,0] = temp1[:,0]*cosz - temp1[:,1]*sinz
temp2[:,1] = temp1[:,0]*sinz + temp1[:,1]*cosz
temp2[:,2] = temp1[:,2]

orbp2[:,0] = temp2[:,0]
orbp2[:,1] = temp2[:,1]*cosx - temp2[:,2]*sinx
orbp2[:,2] = temp2[:,1]*sinx + temp2[:,2]*cosx
 