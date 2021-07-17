# -*- coding: utf-8 -*-
from astropy.time import Time
from astropy.io import fits
from astropy import wcs
import datetime
import easygui
import os
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

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
comveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'eq')
obsveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'obs,eq')
comobs = orb_obs_new(comdenom, obsloc, pysav, orbitdir)

#%%
times = Time(obsveceq[:,0],format='jd')

dir_list = sorted(os.listdir(imagedir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

res = np.zeros(fits_total)
ts = np.zeros(fits_total)

for image_id in range(fits_total):

    #fits name and save file name
    fitstemp = os.path.join(imagedir, fits_list[image_id])
    
    #get fits data
    fitsimg = fits.open(fitstemp)
    try:
        w = wcs.WCS(fitsimg[0].header, key = 'A')
    except:
        w = wcs.WCS(fitsimg[0].header)
    data = fitsimg[0].data
    
    #get image time from filename
    cmin = int(os.path.basename(fitstemp)[11:13])
    chour = int(os.path.basename(fitstemp)[9:11])
    cday = int(os.path.basename(fitstemp)[6:8])
    cmonth = int(os.path.basename(fitstemp)[4:6])
    cyear = int(os.path.basename(fitstemp)[0:4])
    csec = 0
    t_img = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                            chour , cmin, csec))
    
    
    #cell of imagetime in ephs
    imgcell = abs(comveceq[:,0] - t_img.jd).argmin()
    
    #get radec
    radec = pos2radec((comveceq[imgcell,6:9] - obsveceq[imgcell,6:9]),False)
    radec = np.array(radec)
    radec = np.resize(radec,(2,1)).T
    
    #get pixel centre
    pixlocs = w.wcs_world2pix(radec,0)
    pnext = pixlocs + np.array([1,1])
    
    nuradecs = w.wcs_pix2world(pnext,0)
    if abs(nuradecs[0,0]-radec[0,0]) > 180:
        nuradecs = nuradecs + ([360,0])
    
    dtheta = np.linalg.norm(nuradecs - radec)/np.sqrt(2)
    rat = np.linalg.norm(comveceq[imgcell,6:9] - obsveceq[imgcell,6:9])
    
    res[image_id] = np.pi*rat*dtheta/180*1.5e8
    ts[image_id] = t_img.jd
    
#%%
tsj  = Time(ts,format='jd') 

r = np.linalg.norm(comveceq[:,6:9] - obsveceq[:,6:9],axis=1)

plt.plot(tsj.jd,res/40000)
plt.plot(times.jd,r)