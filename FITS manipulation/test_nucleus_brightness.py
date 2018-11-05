#*************************************************************
#Plots a graph of nucleus brightness over a whole dataset
#*************************************************************

import easygui
import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
import datetime
import matplotlib.path as mplPath
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from imagetime_methods import image_time_yudish, image_time_user, image_time_stereo, image_time_filename_yuds
from orbitdata_loading_functions import orb_vector, orb_obs
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
com_in_image = True
fidx = 0

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
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,eq')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'eq,d10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

dir_list = sorted(os.listdir(imagedir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)
   
#%%**********************
#SECOND CELL - get radecs
#************************
nucvals = np.zeros((300),dtype=float)
jdeez = np.zeros((300),dtype=float)

while((fidx < fits_total) and com_in_image == True):
       
    fitsinfile = fits_list[fidx]
    filebase = fits_list[fidx].split('.')[0]
    fitsin = os.path.join(imagedir,fitsinfile)
    
    #ensures image inputted correctly depending on size of data cube
    [colr, colg, colb, fitscoords] = correct_for_imagetype(imagedir, fitsin, fitsinfile)
    
    #get RA/DEC data    
    onedimg = fits.open(fitscoords)
    if 'Stereo' in obsloc:
        if comdenom == 'c2011l4':
            w = wcs.WCS(onedimg[0].header)
        elif comdenom == 'c2006p1':
            if 'A' in obsloc:
                try:
                    w = wcs.WCS(onedimg[0].header, key = 'A')
                except:
                    w = wcs.WCS(onedimg[0].header)   
            else:
                w = wcs.WCS(onedimg[0].header)
    elif 'Soho' in obsloc:
        w = wcs.WCS(onedimg[0].header)
    elif 'Earth' or 'ISS' in obsloc:
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
    
    [ctime,uncertainty_range_exists] = image_time_stereo(filebase)

    #find relevant observer parameters of comet at observer time
    comcel = np.where(abs(obsveceq[:,0] - ctime.jd)==abs(obsveceq[:,0] - ctime.jd).min())[0][0]
    
    LT_cor = int(np.round(np.linalg.norm(comveceq[comcel,6:9] - 
    obsveceq[comcel,6:9])*8.316746397269274))
    com_ra_dec = pos2radec(comveceq[comcel - LT_cor,6:9] - obsveceq[comcel,6:9],fixwrapsbool)
    
    com_in_image = com_Path.contains_point((com_ra_dec[0]-np.mean(ra_m-ra),com_ra_dec[1]))
    
    if com_in_image == False:
        break
    
    com_ra_dec_array = np.zeros((1,2))
    com_ra_dec_array[0,0] = com_ra_dec[0] - np.mean(ra_m-ra); com_ra_dec_array[0,1] = com_ra_dec[1]
    
    compos = tuple(np.round(w.wcs_world2pix(com_ra_dec_array,0)[0][::-1]).astype(int))
    nucvals[fidx] = colr[compos]
    jdeez[fidx] = ctime.jd
    fidx +=1

jdeez = jdeez[0:fidx]
nucvals = nucvals[0:fidx]

plt.plot(jdeez,nucvals)

npsave = os.path.join(os.path.dirname(imagedir), inst + '_thresholds.npy')
np.save(npsave,nucvals)
