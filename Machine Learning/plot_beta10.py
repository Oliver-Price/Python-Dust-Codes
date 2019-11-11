#*************************************************************
#Program to find the time at which an image was taken, with hard coding only
#*************************************************************

import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
import matplotlib.path as mplPath
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import astropy.time
import pickle

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from orbitdata_loading_functions import orb_vector_new, orb_obs_new
from FP_diagnostics import plot_orbit, plot_sunearth_vec
from FP_plot_functions import ra2xpix, dec2ypix, setaxisup, plotpixel
from conversion_routines import pos2radec, fixwraps
from io_methods import correct_for_imagetype
from particle_sim import part_sim

#%%**********************
#FIRST CELL - GET DATA IN
#************************

denoms = ('2P','45P','109P','153P','c1995o1','c1996b2','c1999s4','c2001q4','c2002c1','c2002t7','c2004f4','c2006m4','c2009p1','c2009r1','c2012k1','c2012s1','c2012x1','c2013r1','c2013v5','c2014q2','c2012f6')
cnames = ('Encke','HMP','Swift-Tuttle','Ikeya-Zhang','Hale-Bopp','Hyakutake','LINEAR','NEAT','Ikeya-Zhang','LINEAR','Bradfield','SWAN','Garradd','McNaught','PAN-Starrs','ISON','LINEAR','Lovejoy','Oukaimeden','Lovejoy','Lemmon')

d2c = dict(zip(denoms,cnames))

fitsdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\fits'
nucdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\astro_nucs\good'
baseimgdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\astro_pngs'
orbitdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\orbit_data'
savedir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\beta5'
pysav = r'C:\PhD\Python\Save_Data\machine'
obsloc = 'Earth'

rasave = os.path.join(pysav,"ra.pickle")
decsave = os.path.join(pysav,"dec.pickle")
goodsave = os.path.join(pysav,"goods.pickle")
denomsave = os.path.join(pysav,"denoms.pickle")
      
with open(rasave, 'rb') as handle:
    ras = pickle.load(handle)
    
with open(decsave, 'rb') as handle:
    decs = pickle.load(handle)

with open(goodsave, 'rb') as handle:
    nlocs = pickle.load(handle)

with open(denomsave, 'rb') as handle:
    dnms = pickle.load(handle)

dir_list = sorted(os.listdir(nucdir))
fits_list = [(s.split('.')[0][:-4] + '.fits') for s in dir_list]

fits_total = len(fits_list)

#%%1.5 CELL - start loop

for fid in range(0,fits_total):
    
    print(fid)
    
    fitsinfile = fits_list[fid]
    filebase = fitsinfile[:fitsinfile.find('.')]
    fitsin = os.path.join(fitsdir,fitsinfile)
    
    imgload = os.path.join(baseimgdir, filebase + '_ast.png')
    
    image_save = os.path.join(savedir,filebase + '.png')

    if os.path.exists(image_save) == False:
            
        fit_parts = fitsinfile.split('_')

        cdenom = fit_parts[0]
        if cdenom[0] == 'c':
            cdenom = cdenom.lower()
        else:
            cdenom = cdenom.upper()
            
        #cdenom for naming, comdenom for orbit lookup
    
        comname = d2c[cdenom]
        
        comdenom = dnms[fid]
        
        obsveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'obs,eq')
        comveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'eq')
        comveceq10 = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'eq,d10')
        comobs = orb_obs_new(comdenom, obsloc, pysav, orbitdir)
        
        #ensures image inputted correctly depending on size of data cube
        [colr, colg, colb, fitscoords] = correct_for_imagetype(fitsdir, fitsin, fitsinfile)
        
        #%%**********************
        #SECOND CELL - Reload Image
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
        #Third Cell - sort out Imgtimehdr
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
        
        comcel = nlocs[fid]
        comcel10 = abs(comveceq10[:,0] - obsveceq[comcel,0]).argmin()
        
        c10t = astropy.time.Time(comveceq10[comcel10,0], format = 'jd')
        c1t = astropy.time.Time(comveceq[comcel,0], format = 'jd')
        
        dt = c1t - c10t
        dtmin = int(round(dt.sec/60))
        
        #find ra and dec of comet
        LT_cor = int(np.round(np.linalg.norm(comveceq[comcel,6:9] - 
        obsveceq[comcel,6:9])*8.316746397269274))
        com_ra_dec = pos2radec(comveceq[comcel - LT_cor,6:9] - obsveceq[comcel,6:9],fixwrapsbool)
        
        #check if comet is within image
        com_box_path = np.zeros((2*sum(np.shape(ra_m))-3,2),dtype =float)
        len1 = np.shape(ra_m)[0]; len0 = np.shape(ra_m)[1]
        com_box_path[0:len0,0] = ra_m[0,:]
        com_box_path[0:len0,1] = dec[0,:]
        com_box_path[len0:len0+len1-1,0] = ra_m[1:,-1]
        com_box_path[len0:len0+len1-1,1] = dec[1:,-1]
        com_box_path[len0+len1-1:2*len0+len1-2,0] = ra_m[-1,-2::-1]
        com_box_path[len0+len1-1:2*len0+len1-2,1] = dec[-1,-2::-1]
        com_box_path[2*len0+len1-2:2*len0+2*len1-3,0] = ra_m[-2::-1,0]
        com_box_path[2*len0+len1-2:2*len0+2*len1-3,1] = dec[-2::-1,0]
       		 
        com_Path = mplPath.Path(com_box_path)    
        
        com_in_image = com_Path.contains_point((com_ra_dec[0],com_ra_dec[1]))
        
        #account for sometimes negative axisdata
        if (axisdata[6] < 0):
            ra_img_lower = axisdata[6] + 360
        else: ra_img_lower = axisdata[6]
        if (axisdata[7] < 0):
            ra_img_higher = axisdata[7] + 360
        else: ra_img_higher = axisdata[7]
        
        backgr_fill = (0,0,0,255)
        featur_fill = (255,255,255,255)
        trajfill = (0,255,255,255)
        trajucfill = (0,0,255,255)
        
        fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
        fnt = ImageFont.truetype(fontloc, 20)  
          
        [ltcomcel, vtraj, vtrajcel] = plot_orbit(comobs,comveceq,obsveceq,axisdata,d,comcel,trajfill,ra_img_lower,
        ra_img_higher,border,pixwidth,rafmin,scale,pixheight,decmin,fnt,featur_fill,com_in_image,fixwrapsbool)      
    
        comsunfill = (0,255,0,255)
    
        plot_sunearth_vec(d,comveceq,obsveceq,ltcomcel,axisdata,ra_img_higher,ra_img_lower,
                    border,pixwidth,pixheight,rafmin,decmin,scale,comsunfill,featur_fill,fnt,fixwrapsbool)
        
        tno = 100
        tvals = np.linspace(0.1,100,tno)
        tidx = 0
        simres = np.zeros((tno,14),dtype = float)
        efinp = obsveceq[comcel,6:9]
        
        while (tidx < tno):
            point_in_image = 1
            simt10min = int(round(144*tvals[tidx]))
            pstart = comveceq10[comcel10-simt10min,6:12]
            sim = part_sim(5,simt10min,10,1,pstart,efinp,dtmin)
            simres[tidx,0] = float(simt10min)/144
            simres[tidx,1] = sim[0] #length of simulation in minutes
            simres[tidx,2] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
            simres[tidx,3:9] = sim[1] #finishing pos/vel
            simres[tidx,9:11] = pos2radec(sim[1][0:3] - 
                obsveceq[int(simres[tidx,2]),6:9],fixwrapsbool)
            point_in_image = int(com_Path.contains_point(
            (simres[tidx,9],simres[tidx,10])))
            simres[tidx,11] = ra2xpix(simres[tidx,9],border,pixwidth,rafmin,scale)
            simres[tidx,12] = dec2ypix(simres[tidx,10],border,pixheight,decmin,scale)                 
            simres[tidx,13] = point_in_image
            tidx += 1
    
        chrfill = (255,192,0,255)
        for ta in range(int(np.sum(simres[:,13]))-1):
            d.line([(simres[ta,11],simres[ta,12]), \
            (simres[ta+1,11],simres[ta+1,12])],\
            fill = chrfill)
        
        comimg.save(image_save)