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
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag)

#%%
        
#z = -87.41489440044852
#x = 77.83726262785106
z = 180 - an
x = inc

cosz = np.cos(np.radians(z))
sinz = np.sin(np.radians(z))
cosx = np.cos(np.radians(x))
sinx = np.sin(np.radians(x))

temp = np.empty_like(comveceq[:,6:9])
orbp = np.empty_like(comveceq[:,6:9])

temp[:,0] = comveceq[:,6]*cosz - comveceq[:,7]*sinz
temp[:,1] = comveceq[:,6]*sinz + comveceq[:,7]*cosz
temp[:,2] = comveceq[:,8]

orbp[:,0] = temp[:,0]
orbp[:,1] = temp[:,1]*cosx - temp[:,2]*sinx
orbp[:,2] = temp[:,1]*sinx + temp[:,2]*cosx
 
        
#%%

cross_vec = np.cross(comveceq[0,6:9],comveceq[-1,6:9])
pn = cross_vec/np.linalg.norm(cross_vec)

#%%
tempv = np.empty_like(comveceq[:,6:9])
orbv = np.empty_like(comveceq[:,6:9])

tempv[:,0] = comveceq[:,9]*cosz - comveceq[:,10]*sinz
tempv[:,1] = comveceq[:,9]*sinz + comveceq[:,10]*cosz
tempv[:,2] = comveceq[:,11]

orbv[:,0] = tempv[:,0]
orbv[:,1] = tempv[:,1]*cosx - tempv[:,2]*sinx
orbv[:,2] = tempv[:,1]*sinx + tempv[:,2]*cosx

#%%

theta_hat = np.empty_like(comveceq[:,6:8])

rho = np.sqrt(orbp[:,0]**2 + orbp[:,1]**2)
theta_hat[:,0] = -orbp[:,1]/rho[:]
theta_hat[:,1] = orbp[:,0]/rho[:]

#%%

v_theta = orbv[:,0]*theta_hat[:,0] + orbv[:,1]*theta_hat[:,1]