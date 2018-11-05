#*************************************************************
#Program to find the time at which an image was taken, with hard coding only
#*************************************************************

import easygui
import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
import datetime
import matplotlib.path as mplPath
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from orbitdata_loading_functions import orb_obs_new
from FP_plot_functions import ra2xpix, dec2ypix, setaxisup, plotpixel
from conversion_routines import pos2radec, fixwraps, find_largest_nonzero_block 
from io_methods import correct_for_imagetype, get_obs_loc, get_hih_low 
from io_methods import get_stereo_instrument, get_soho_instrument

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
    timedir = cdata[28][26:-2]
    pysav = cdata[27][24:-2]
    horiztag = cdata[40][10:]
    obslocstr = cdata[34][19:]
#choose observer locations
[obsloc, imagedir] = get_obs_loc(obslocstr, imagedir)
if "Stereo" in obsloc: [inst, imagedir] = get_stereo_instrument(imagedir)
elif obsloc == "Soho": [inst, imagedir] = get_soho_instrument(imagedir)
else: inst = ''  
#import the orbit data
comobs = orb_obs_new(comdenom, obsloc, pysav, orbitdir, horiztag)

#choosing fits file to display and getting pathnames
fitsin = easygui.fileopenbox(default = os.path.join(imagedir,'*'))
if fitsin == '.': sys.exit("No file selected")
fitsinfile = os.path.basename(fitsin)
filebase = fitsinfile[:fitsinfile.find('.')]

#ensures image inputted correctly depending on size of data cube
[colr, colg, colb, fitscoords] = correct_for_imagetype(imagedir, fitsin, fitsinfile)

        
#%%**********************
#SECOND CELL - Plot Image
#************************

#get RA/DEC data    
onedimg = fits.open(fitscoords)
if 'Stereo' in obsloc:
    if comdenom == 'c2011l4':
        w = wcs.WCS(onedimg[0].header)
    elif comdenom == 'c2006p1':
        w = wcs.WCS(onedimg[0].header)#, key = 'A')
    if 'diff' in inst or 'MGN' in inst: plotmethodlog = False
    else: plotmethodlog = True
elif 'Soho' in obsloc:
    w = wcs.WCS(onedimg[0].header)
    if 'diff' in inst or 'MGN' in inst: plotmethodlog = False
    else: plotmethodlog = True
elif 'Earth' or 'ISS' in obsloc:
    plotmethodlog = False
    w = wcs.WCS(onedimg[0].header)

#make a 2xN array of all pixel locations
ya = onedimg[0].data.shape[0]
xa = onedimg[0].data.shape[1]
xv, yv = np.meshgrid(np.arange(xa), np.arange(ya))
xv = np.reshape(xv, (xa*ya,1))
yv = np.reshape(yv, (xa*ya,1))
coords = np.zeros((xa*ya,2))
coords[:,0] = xv[:,0]
coords[:,1] = yv[:,0]

#convert each pixel locaion to an RA and DEC array
radecs = w.wcs_pix2world(coords, 0)
ra = np.reshape(radecs[:,0], (ya,xa))
dec = np.reshape(radecs[:,1], (ya,xa))
    
#find minimum/maximum values
ramin = np.amin(ra)
ramax = np.amax(ra)
decmin = np.amin(dec)
decmax = np.amax(dec)
    
[ra_m, rafmin, rafmax, fixwrapsbool] = fixwraps(ra, ramax, ramin)
    
#make a canvas with a fixed pixel height and border
pixheight = 800
pixwidth = int(pixheight*(rafmax - rafmin)/(decmax - decmin))
border = 100

scale = pixheight/(decmax - decmin)
imgwidth = pixwidth+int(4*border)
imgheight = pixheight+int(3*border)

#%%**************************
#FOURTH CELL - Get Imgtimehdr
#****************************

#initialise comet box path
com_box_path = np.zeros((2*sum(np.shape(ra))-3,2),dtype =float)
len1 = np.shape(ra)[0]; len0 = np.shape(ra)[1]
com_box_path[0:len0,0] = ra[0,:]
com_box_path[0:len0,1] = dec[0,:]
com_box_path[len0:len0+len1-1,0] = ra[1:,-1]
com_box_path[len0:len0+len1-1,1] = dec[1:,-1]
com_box_path[len0+len1-1:2*len0+len1-2,0] = ra[-1,-2::-1]
com_box_path[len0+len1-1:2*len0+len1-2,1] = dec[-1,-2::-1]
com_box_path[2*len0+len1-2:2*len0+2*len1-3,0] = ra[-2::-1,0]
com_box_path[2*len0+len1-2:2*len0+2*len1-3,1] = dec[-2::-1,0]
com_Path = mplPath.Path(com_box_path)   

comlocra = ra2xpix(comobs[:,5], border, pixwidth, rafmin, scale)
comlocdec = dec2ypix(comobs[:,6], border, pixheight, decmin, scale)

decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)

traj_pix = w.wcs_world2pix(comobs[goodlocs,5:7],0)
locs_1 = np.where((traj_pix[:,1]>0)&(traj_pix[:,0]>0))
traj_pix = traj_pix[locs_1[0],:]
locs_2 = np.where((traj_pix[:,0]<(colr.shape[0]-1))&(traj_pix[:,1]<(colr.shape[1]-1)))
traj_pix = traj_pix[locs_2[0],:]

traj_pix_floor = np.floor(traj_pix).astype(int)
traj_pix_round = np.round(traj_pix).astype(int)

trajvals = np.zeros((traj_pix[:,0].size,4))

for p in range(locs_2[0].size):
    
    x = traj_pix[p,0]
    y = traj_pix[p,1]
    x_r = traj_pix[p,0]
    y_r = traj_pix_round[p,1]
    x_f = traj_pix_floor[p,0]
    y_f = traj_pix_floor[p,1]
    trajvals[p,0] = (colr[y_f,x_f]*(y_f - y + 1)*(x_f - x + 1) + 
                   colr[y_f+1,x_f]*(y - y_f)*(x_f - x + 1) +
                   colr[y_f,x_f+1]*(y_f - y + 1)*(x - x_f) + 
                   colr[y_f+1,x_f+1]*(y - y_f)*(x - x_f))
    trajvals[p,1] = (colg[y_f,x_f]*(y_f - y + 1)*(x_f - x + 1) + 
                   colg[y_f+1,x_f]*(y - y_f)*(x_f - x + 1) +
                   colg[y_f,x_f+1]*(y_f - y + 1)*(x - x_f) + 
                   colg[y_f+1,x_f+1]*(y - y_f)*(x - x_f))
    trajvals[p,2] = (colb[y_f,x_f]*(y_f - y + 1)*(x_f - x + 1) + 
                   colb[y_f+1,x_f]*(y - y_f)*(x_f - x + 1) +
                   colb[y_f,x_f+1]*(y_f - y + 1)*(x - x_f) + 
                   colb[y_f+1,x_f+1]*(y - y_f)*(x - x_f))
    trajvals[p,3] = np.mean(trajvals[p,0:3])

