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
from scipy.signal import find_peaks

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from orbitdata_loading_functions import orb_obs_new
from FP_plot_functions import ra2xpix, dec2ypix, setaxisup, plotpixel
from conversion_routines import pos2radec, fixwraps
from io_methods import correct_for_imagetype, get_obs_loc, get_hih_low 
from io_methods import get_stereo_instrument, get_soho_instrument

#%%**********************
#FIRST CELL - GET DATA IN
#************************

denoms = ('2P','45P','109P','153P','c1995o1','c1996b2','c1999s4','c2001q4','c2002c1','c2002t7','c2004f4','c2006m4','c2009p1','c2009r1','c2012k1','c2012s1','c2012x1','c2013r1','c2013v5','c2014q2','c2012f6')
cnames = ('Encke','HMP','Swift-Tuttle','Ikeya-Zhang','Hale-Bopp','Hyakutake','LINEAR','NEAT','Ikeya-Zhang','LINEAR','Bradfield','SWAN','Garradd','McNaught','PAN-Starrs','ISON','LINEAR','Lovejoy','Oukaimeden','Lovejoy','Lemmon')

d2c = dict(zip(denoms,cnames))

fitsdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\fits'
imgdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\astro_pngs'
nucdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\astro_nucs'
orbitdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\orbit_data'
pltsav = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\Orbitpathplots60+200'
pysav = r'C:\PhD\Python\Save_Data\machine'
obsloc = 'Earth'

bads = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\astro_nucs\bad'

good = 0
bad = 0

if not os.path.exists(nucdir): os.makedirs(nucdir)

dir_list = sorted(os.listdir(imgdir))
fits_list = [(s.split('.')[0][:-4] + '.fits') for s in dir_list]

fits_total = len(fits_list)

#%%1.5 CELL - start loop

for fid in range(0,fits_total):
    
    print(fid)
    
    fitsinfile = fits_list[fid]
    filebase = fitsinfile[:fitsinfile.find('.')]
    fitsin = os.path.join(fitsdir,fitsinfile)
    
    fit_parts = fitsinfile.split('_')
    
    comdenom = fit_parts[0]
    if comdenom[0] == 'c':
        comdenom = comdenom.lower()
    else:
        comdenom = comdenom.upper()

    comname = d2c[comdenom]
    
    #ensures image inputted correctly depending on size of data cube
    [colr, colg, colb, fitscoords] = correct_for_imagetype(fitsdir, fitsin, fitsinfile)
    
    #%%**********************
    #SECOND CELL - Plot Image
    #************************
    
    #get RA/DEC data    
    onedimg = fits.open(fitscoords)
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
    com_box_path[0:len0,0] = ra_m[0,:]
    com_box_path[0:len0,1] = dec[0,:]
    com_box_path[len0:len0+len1-1,0] = ra_m[1:,-1]
    com_box_path[len0:len0+len1-1,1] = dec[1:,-1]
    com_box_path[len0+len1-1:2*len0+len1-2,0] = ra_m[-1,-2::-1]
    com_box_path[len0+len1-1:2*len0+len1-2,1] = dec[-1,-2::-1]
    com_box_path[2*len0+len1-2:2*len0+2*len1-3,0] = ra_m[-2::-1,0]
    com_box_path[2*len0+len1-2:2*len0+2*len1-3,1] = dec[-2::-1,0]
    com_Path = mplPath.Path(com_box_path)   

    #import the orbit data
    try:
        comobs = orb_obs_new(comdenom, obsloc, pysav, orbitdir)
    except:
        comobs = orb_obs_new(comdenom + '_1', obsloc, pysav, orbitdir)
        
    cra = comobs[:,5]
    
    if cra.ptp() > 359:
        cra[cra<180] = cra[cra<180] + 360
    
    comlocra = ra2xpix(cra, border, pixwidth, rafmin, scale)
    comlocdec = dec2ypix(comobs[:,6], border, pixheight, decmin, scale)
    
    decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
    ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
    goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)
    
    if np.size(goodlocs) == 0:
        try:
            comobs = orb_obs_new(comdenom + '_2', obsloc, pysav, orbitdir)
            
            cra = comobs[:,5]
            
            if cra.ptp() > 359:
                cra[cra<180] = cra[cra<180] + 360
            
            comlocra = ra2xpix(cra, border, pixwidth, rafmin, scale)
            comlocdec = dec2ypix(comobs[:,6], border, pixheight, decmin, scale)
            
            decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
            ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
            goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)
        except:
            pass
        
    if np.size(goodlocs) == 0:
        try:
            comobs = orb_obs_new(comdenom + '_3', obsloc, pysav, orbitdir)
            
            cra = comobs[:,5]
            
            if cra.ptp() > 359:
                cra[cra<180] = cra[cra<180] + 360
            
            comlocra = ra2xpix(cra, border, pixwidth, rafmin, scale)
            comlocdec = dec2ypix(comobs[:,6], border, pixheight, decmin, scale)
            
            decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
            ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
            goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)
        except:
            pass
        
    #%%**************************
    #FIFTH CELL - Get relavent traj points
    #****************************

    try:
        goodradecs = comobs[goodlocs,5:7]
        goodradecs_t = tuple(map(tuple, goodradecs))
        goodbois =  np.array(com_Path.contains_points(goodradecs_t))
        goodlocs = goodlocs[np.nonzero(goodbois*goodlocs)]
        
        traj_pix = w.wcs_world2pix(comobs[goodlocs,5:7],0)
        traj_pix_floor = np.floor(traj_pix).astype(int)
        traj_pix_round = np.round(traj_pix).astype(int)
        
        trajvals = np.zeros((traj_pix[:,0].size,4))
        
        for p in range(goodlocs.size):
            
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
        
        #%%**************************
        #SIXTH CELL - Plot stuff
        #****************************
        wu = 20
        x = trajvals[:,3]
        peaks, _ = find_peaks(x, width = int(goodlocs.size/wu), height = 180)
        
        while (peaks.size == 0) and (wu < 300):
            wu += 20
            x = trajvals[:,3]
            peaks, _ = find_peaks(x, width = int(goodlocs.size/wu), height = 200)
        
        h = 200
        while (peaks.size == 0) and (h > 150):
            wu = 200
            h -= 10
            x = trajvals[:,3]
            peaks, _ = find_peaks(x, width = int(goodlocs.size/wu), height = h)
            
        if peaks.size > 1:
            peaks = np.array([int(round(peaks.mean()))])
            
        iload = os.path.join(imgdir,filebase + '_ast.png')

        comimg = Image.open(iload)
        d = ImageDraw.Draw(comimg)
    
        xl = ra2xpix(cra[goodlocs][peaks[0]], border, pixwidth, rafmin, scale)
        yl = dec2ypix(comobs[goodlocs,6][peaks[0]], border, pixheight, decmin, scale)
    
        xsiz = 5
        d.line( [ ( xl - xsiz , yl - xsiz ) ,
          ( xl + xsiz , yl + xsiz ) ] ,
          fill = (255,0,255,255) )  
        d.line( [ ( xl - xsiz , yl + xsiz ) ,
          ( xl + xsiz , yl - xsiz ) ] ,
          fill = (255,0,255,255) )
        
        isave = os.path.join(nucdir,filebase+ '_nuc.png')
        
        comimg.save(isave)
        comimg.close()
            
    except:
        pass
    
print(good,bad)