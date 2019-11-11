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
from astropy.time import Time
import matplotlib.pyplot as plt

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from orbitdata_loading_functions import *
from io_methods import *

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

#import the orbit data
com_vec = orb_vector(comdenom, 'Stereo_A', pysav, orbitdir,
                      horiztag)
sterb_vec = orb_vector(comdenom, 'Soho', pysav, orbitdir,
                      horiztag, opts = 'obs')
earth_vec = orb_vector(comdenom, 'Earth', pysav, orbitdir,
                      horiztag, opts = 'obs')

#%%
x = 38989
f, ax = plt.subplots(1)
#ax.scatter(stera_vec[x,6],stera_vec[x,7],c='r',marker='X')
#ax.scatter(sterb_vec[x,6],sterb_vec[x,7],c='g',marker='X')
ax.scatter(earth_vec[x,6],earth_vec[x,7],c='b',marker='X')
ax.scatter(0,0,c='k',marker='X')
ax.scatter(com_vec[x,6],com_vec[x,7],c='m',marker='X')
ax.set_xlim([-1.5,1.5])
ax.set_ylim([-1.5,1.5])
ax.set_aspect('equal')
