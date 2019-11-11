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
from conversion_routines import pos2radec, fixwraps
from io_methods import correct_for_imagetype, get_obs_loc, get_hih_low 
from io_methods import get_stereo_instrument, get_soho_instrument

#%%**********************
#FIRST CELL - GET DATA IN
#************************

denoms = ('2P','45P','109P','153P','c1995o1','c1996b2','c1999s4','c2001q4','c2002c1','c2002t7','c2004f4','c2006m4','c2009p1','c2009r1','c2012k1','c2012s1','c2012x1','c2013r1','c2013v5','c2014q2','c2012f6')
cnames = ('Encke','HMP','Swift-Tuttle','Ikeya-Zhang','Hale-Bopp','Hyakutake','LINEAR','NEAT','Ikeya-Zhang','LINEAR','Bradfield','SWAN','Garradd','McNaught','PAN-Starrs','ISON','LINEAR','Lovejoy','Oukaimeden','Lovejoy','Lemmon')

d2c = dict(zip(denoms,cnames))

imagedir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\fits'
imgsavedir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\astro_pngs'
orbitdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\orbit_data'
pysav = r'C:\PhD\Python\Save_Data\machine'
obsloc = 'Earth'

dir_list = sorted(os.listdir(imagedir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

#%%1.5 CELL - start loop

for fid in range(22,23):#fits_total):
    
    fitsinfile = fits_list[fid]
    filebase = fitsinfile[:fitsinfile.find('.')]
    fitsin = os.path.join(imagedir,fitsinfile)
    
    imgsave = os.path.join(imgsavedir, filebase + '_ast.png')

    if True:#os.path.exists(imgsave) == False:
            
        fit_parts = fitsinfile.split('_')
        
        comdenom = fit_parts[0]
        if comdenom[0] == 'c':
            comdenom = comdenom.lower()
        else:
            comdenom = comdenom.upper()
    
        comname = d2c[comdenom]
        
        #ensures image inputted correctly depending on size of data cube
        [colr, colg, colb, fitscoords] = correct_for_imagetype(imagedir, fitsin, fitsinfile)
    
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
        
        hih = 255
        low = 0
        
        scale = pixheight/(decmax - decmin)
        imgwidth = pixwidth+int(4*border)
        imgheight = pixheight+int(3*border)
        comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,(0,0,0,255))
        d = ImageDraw.Draw(comimg)
        
        colcr = np.clip(255.0/(hih-low)*(colr-low),0,255).astype(int)
        colcg = np.clip(255.0/(hih-low)*(colg-low),0,255).astype(int)
        colcb = np.clip(255.0/(hih-low)*(colb-low),0,255).astype(int)         
        
        for x in range(0,np.shape(colr)[0]-1):
            for y in range (0,np.shape(colr)[1]-1):
                plotpixel(d,x,y,ra_m,dec,border,pixwidth,pixheight,decmin,rafmin,
                          scale,colcr[x,y],colcg[x,y],colcb[x,y])
                    
        #%%********************
        #THIRD CELL - Draw Axis
        #**********************
        featur_fill = (255,255,255,255)
        #draws a border       
        a = d.polygon([(border,border),(border*2+pixwidth,border), \
        (border*2+pixwidth,border*2+pixheight),(border,border*2+pixheight)], \
        outline = featur_fill)
        
        #most of the dirty stuff is bunged into this function
        axisdata = setaxisup(rafmax,rafmin,decmax,decmin,border,pixheight,pixwidth,scale)
        
        majt = 20  #major tick length
        mint = 10  #minor tick length
        rdaxis = pixheight + border*2
        
        fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
        fnt = ImageFont.truetype(fontloc, 20)
        smallfnt = ImageFont.truetype(fontloc, 10)
        largefnt = ImageFont.truetype(fontloc, 30)
        
        for div in range(0, (np.size(axisdata[1]))): #RA axis major ticks
            b = d.line([(axisdata[1][div],rdaxis-majt),(axisdata[1][div],rdaxis)],\
            fill = featur_fill)
            tick = "{:.1f}".format(axisdata[0][div]%360)
            d.text((axisdata[1][div] - len(tick)*5,rdaxis + 10), \
            tick, font=fnt, fill= featur_fill)
            
        for div in range(0, (np.size(axisdata[2]))): #RA axis minor ticks
            b = d.line([(axisdata[2][div],rdaxis-mint),(axisdata[2][div],rdaxis)],\
            fill= featur_fill)
        
        for div in range(0, (np.size(axisdata[4]))): #DEC axis major ticks
            b = d.line([(border+majt,axisdata[4][div]),(border,axisdata[4][div])],\
            fill= featur_fill)
            tick = "{:.1f}".format(axisdata[3][div])
            d.text((border - len(tick)*5 - 40,axisdata[4][div] - 10 ), \
            tick, font=fnt, fill=(255,255,255,128))
            
        for div in range(0, (np.size(axisdata[5]))): #DEC axis minor ticks
            b = d.line([(border+mint,axisdata[5][div]),(border,axisdata[5][div])],\
            fill= featur_fill)
        
        #axis labels
        d.text((1.5*border + pixwidth*0.5 - 145,pixheight + 2.5*border), \
        "Right Ascension (Degrees)", font=fnt, fill= featur_fill)
        d.text((0.25*border - 10,0.75*border - 20), \
        "Declination (Degrees)", font=fnt, fill= featur_fill)
        
        #plot title
        plttitle = (comdenom.upper() + ' ' + comname)
        d.text((1.5*border + pixwidth*0.5 - len(plttitle)*5 - 60,.35*border), \
        plttitle, font=largefnt, fill= featur_fill)
        
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
        
        for o in range(0,np.size(goodlocs)-1):
            ocel = goodlocs[o]
            d.line( [ (comlocra[ocel] ,comlocdec[ocel]), 
            (comlocra[ocel+1] ,comlocdec[ocel+1]) ],
            fill = (255,0,0,255))
    
        comimg.save(imgsave,'png')
    
        print ("Done!")