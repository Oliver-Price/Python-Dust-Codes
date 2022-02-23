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
import datetime

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

#choose observer locations
[obsloc, imagedir] = get_obs_loc(obslocstr, imagedir)
if "Stereo" in obsloc: [inst, imagedir] = get_stereo_instrument(imagedir)
elif obsloc == "Soho": [inst, imagedir] = get_soho_instrument(imagedir)
else: inst = ''
 
#import the orbit data
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag)
#comobs = orb_obs_new(comdenom, obsloc, pysav, orbitdir)

#%%
vr = np.linalg.norm(comveceq[:,9:12],axis=1)*1.731e+3
print(vr.max())

#%%
times = Time(obsveceq[:,0],format='jd')

t_img1 = astropy.time.Time(datetime.datetime(2013, 3, 12, 5, 0, 0))
t_img2 = astropy.time.Time(datetime.datetime(2013, 3, 13, 17, 30, 0))

tcell1 = abs(comveceq[:,0] - t_img1.jd).argmin()
tcell2 = abs(comveceq[:,0] - t_img2.jd).argmin()

r = np.linalg.norm(comveceq[:,6:9],axis=1)
rs = np.linalg.norm(comveceq[:,6:9] - obsveceq[:,6:9],axis=1)