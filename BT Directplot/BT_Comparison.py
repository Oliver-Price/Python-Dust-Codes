# -*- coding: utf-8 -*-
#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************

import string
import easygui
import os
import sys
import pickle
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont
import webbrowser

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")  

from orbitdata_loading_functions import orb_vector, orb_obs
from BT_plot_functions import beta2ypix, simt2xpix, plotpixel, logbeta2ypix, logplotpixel, bt_setaxisup
from io_methods import get_obs_loc, get_stereo_instrument, get_soho_instrument, get_hih_low
from conversion_routines import fixwraps, round_to_base

#%%**********************
#FIRST CELL - GET DATA IN
#************************

#choosing comet data to use
inputfilefolder = "C:\PhD\Comet_data\Input_files\*pt1.txt"
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
                      horiztag, opts = 'obs,eq')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'eq,d10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)
