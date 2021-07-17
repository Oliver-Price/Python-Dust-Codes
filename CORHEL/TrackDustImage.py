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

#%%

#dir properties
fitsdir = r'C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\HI-1-diff'
pngdir = os.path.join(fitsdir,'dust_tracked2')
if not os.path.exists(pngdir): os.makedirs(pngdir)

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)


#draw properties
low = -1800; hih = 1700
methodlog = False

#dust properties
t_dust = Time('2013-03-05T12:00:00.000',format='isot')
beta = 0.8

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
obsveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'obs,eq')

#sim properties
startcell = abs(comveceq[:,0] - t_dust.jd).argmin()
simt = 20*60*24 #10 days in minutes
nperday = 1440
tstart = comveceq[startcell,6:12]

#simulating
(times,positions) = part_sim_fine_track(beta, simt, nperday, tstart)
jds = Time(t_dust.jd+times/1440,format='jd').squeeze()

poslist = np.empty((6,0))
jdlist = np.empty((1,0))

for fits_id in range(0,fits_total):
        
        #fits name and save file name
        fitstemp = os.path.join(fitsdir, fits_list[fits_id])
        pngtemp = os.path.join(pngdir, (fits_list[fits_id].split(".")[0] + "dust.png"))

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
        endcell = abs(comveceq[:,0] - t_img.jd).argmin()
        
        #position of dust at imagetime
        dcell = abs((jds-t_img).jd).argmin()
        position = positions[dcell]
        
        #diagnostics
        poslist = np.append(poslist,np.expand_dims(position, axis=1),axis = 1)
        jdlist = np.append(jdlist,t_img.jd)
        
        #get radec
        radec = pos2radec((position[0:3] - obsveceq[endcell,6:9]),False)
        radec = np.array(radec)
        radec = np.resize(radec,(2,1)).T
        
        #get fits data
        fitsimg = fits.open(fitstemp)
        w = wcs.WCS(fitsimg[0].header)
        data = fitsimg[0].data
        
        #get pixel centre
        pixlocs = w.wcs_world2pix(radec,0)
        pround = np.floor(pixlocs)[0].astype(int)
        
        #check if pixel centre in image
        imagegood = (0 <= pround[0] <= 1024)&(0 <= pround[1] <= 1024)
        
        if imagegood == True:
            #get 100x100 area to cut
            xstart = pround[1] - 50
            xend = pround[1] + 51
            ystart = pround[0] - 50
            yend = pround[0] + 51
            
            #get start and end of area to cut from image
            #get start and end offset of area to paste in center of new image
            if xstart < 0:
                xoff = -xstart
                xstart = 0
            else:
                xoff = 0
            if xend > data.shape[0]:
                xend = data.shape[0]
            xoff2 = xend-xstart+xoff
                        
            if ystart < 0:
                yoff = -ystart
                ystart = 0
            else:
                yoff = 0
            if yend > data.shape[1]:
                yend = data.shape[1]
            yoff2 = yend-ystart+yoff
               
            #cut dat
            cut = data[xstart:xend,ystart:yend]
            cut = np.clip((255.0/(hih - low)*(cut-low)),0,255).astype(int)
            
            #paste data, flip, resize to 3x and add rgb layers
            pre = np.full((101,101),128)
            pre[xoff:xoff2,yoff:yoff2] = cut
            pre = np.flipud(pre)
            pre = pre.reshape(pre.shape[0],pre.shape[1],1)
            pre = pre.repeat(3,axis=0).repeat(3,axis=1).repeat(3,axis=2)
            
            #draw cross at image centre
            xc = round(pre.shape[0]/2)
            yc = round(pre.shape[1]/2)
            clen = 7
            purp = np.array([0,255,0])
            for xct in [xc-1,xc,xc+1]:
                for yct in [yc-1,yc,yc+1]:
                    for l in range(clen):
                        pre[xct+l,yct+l,:] = purp
                        pre[xct+l,yct-l,:] = purp
                        pre[xct-l,yct-l,:] = purp
                        pre[xct-l,yct+l,:] = purp
              
            #save
            matplotlib.image.imsave(pngtemp, pre.astype('uint8'))