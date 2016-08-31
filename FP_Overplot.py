#*************************************************************
#Program to visualise fits image of comet and overlay a
#finson-probstein diagram according to it's orbital parameters
#*************************************************************

import urllib2
import string
import easygui
import os
import numpy as np
import math as m
from scipy.io.idl import readsav
import sys
from astropy.io import fits
from astropy import wcs
import astropy.time
import datetime
import pickle
import matplotlib.path as mplPath
from PIL import Image, ImageDraw, ImageFont
from orbitdata_loading_functions import orb_vector, orb_obs
from FP_plot_functions import ra2xpix, dec2ypix, setaxisup, plotpixel
from FP_plot_functions import draw_synchrones, draw_sydynes, draw_datap,
from FP_plot_functions import draw_data_reg, annotate_plotting
from FP_diagnostics import plot_orbit, plot_orbit_points, plot_sunearth_vec
from FP_diagnostics import write_bt_ranges, write_properties, plot_compos_unc
from imagetime_methods import image_time_yudish, image_time_user, image_time_stereo
from conversion_routines import pos2radec, fixwraps, correct_for_imagetype, 
from conversion_routines import col_corrections, get_obs_loc
from simulation_setup import simulation_setup
from particle_sim import part_sim

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
    
#import the orbit data
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,eq')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'eq,d10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

#choosing fits file to display and getting pathnames
fitsin = easygui.fileopenbox(default = os.path.join(imagedir, r"*.fits"),
                             filetypes= r"*.fits")
fitsinfile = os.path.basename(fitsin)
filebase = fitsinfile[:string.find(fitsinfile,'.')]

inv = False
#check if image exists already/invert it
if inv == False:   
    imgsav = os.path.join(imagedir, filebase + '.png')
    imgexists = os.path.exists(imgsav)
elif inv == True:
    imgsav = os.path.join(imagedir, filebase + '_inverted.png')
    imgexists = os.path.exists(imgsav)  

#parameter savefile locations
picklesavefile = os.path.join(pysav, filebase + '_dustplot')
picklexists = os.path.exists(picklesavefile)

forceredraw = True
if (imgexists == False) or (picklexists == False) or (forceredraw == True):
    
    #ensures image inputted correctly depending on size of data cube
    [colr, colg, colb, fitscoords] = correct_for_imagetype(imagedir, fitsin, fitsinfile)
    
    #correcting for wrong data type, correct range to 255, account for inversion
    [backgr_fill, featur_fill, colcr, colcg, colcb] = col_corrections(inv,colr,colg,colb)
        
#%%**********************
#SECOND CELL - Plot Image
#************************

    #get RA/DEC data    
    onedimg = fits.open(fitscoords)
    if 'Earth' in obsloc:
        w = wcs.WCS(onedimg[0].header)
    elif 'Stereo' in obsloc:
        w = wcs.WCS(onedimg[0].header, key = 'A')
    elif 'Soho' in obsloc:
        sys.exit("soho data not yet implemented")
    
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
    radecs = w.wcs_pix2world(coords, 1)
    ra = np.reshape(radecs[:,0], (ya,xa))
    dec = np.reshape(radecs[:,1], (ya,xa))
    
    #find minimum/maximum values
    ramin = np.amin(ra)
    ramax = np.amax(ra)
    decmin = np.amin(dec)
    decmax = np.amax(dec)
    
    [ra_m, rafmin, rafmax] = fixwraps(ra, ramax, ramin)
    
    #make a canvas with a fixed pixel height and border
    pixheight = 800
    pixwidth = int(pixheight*(rafmax - rafmin)/(decmax - decmin))
    border = 100
    scale = pixheight/(decmax - decmin)
    imgwidth = pixwidth+int(4*border)
    imgheight = pixheight+int(3*border)
    comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,backgr_fill)
    d = ImageDraw.Draw(comimg)
    
    if obsloc != "Earth":
        plotmethodlog = True
    else: plotmethodlog = False
    
    if plotmethodlog == True:
        if comdenom == 'c2011l4':
            low = 3000
            hih = 20000
        elif comdenom == 'c2006p1':
            low = 10000
            hih = 1500000        
        
        for x in xrange(0, ya-2):
            for y in xrange(0, xa-2):
                fillval = sorted([1, colr[x,y], 9999999999])[1]
                fillco = int(round(255*(np.log10(fillval) - np.log10(low))*
    								1 / (np.log10(hih) - np.log10(low))))
                fillco = sorted([0, fillco, 255])[1]
                plotpixel(d,x,y,ra_m,dec,border,pixwidth,pixheight,decmin,rafmin,
                          scale,fillco,fillco,fillco)
                
    elif plotmethodlog == False:            
        for x in xrange(0, ya-2):
            for y in xrange(0, xa-2):
                plotpixel(d,x,y,ra_m,dec,border,pixwidth,pixheight,decmin,rafmin,
                scale,colcr[x,y],colcg[x,y],colcb[x,y])
            
#%%********************
#THIRD CELL - Draw Axis
#**********************
    
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
    
    for div in xrange(0, (np.size(axisdata[1]))): #RA axis major ticks
        b = d.line([(axisdata[1][div],rdaxis-majt),(axisdata[1][div],rdaxis)],\
        fill = featur_fill)
        tick = str(axisdata[0][div]%360)
        d.text((axisdata[1][div] - len(tick)*5,rdaxis + 10), \
        tick, font=fnt, fill= featur_fill)
        
    for div in xrange(0, (np.size(axisdata[2]))): #RA axis minor ticks
        b = d.line([(axisdata[2][div],rdaxis-mint),(axisdata[2][div],rdaxis)],\
        fill= featur_fill)
    
    for div in xrange(0, (np.size(axisdata[4]))): #DEC axis major ticks
        b = d.line([(border+majt,axisdata[4][div]),(border,axisdata[4][div])],\
        fill= featur_fill)
        tick = str(axisdata[3][div])
        d.text((border - len(tick)*5 - 40,axisdata[4][div] - 10 ), \
        tick, font=fnt, fill=(255,255,255,128))
        
    for div in xrange(0, (np.size(axisdata[5]))): #DEC axis minor ticks
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
    
    pix = comimg.load()
    
#%%**************************
#FOURTH CELL - Get Imgtimehdr
#****************************
    
    #some amateur images lack times, need to sometimes use yudish imgtimehdr
    if obsloc == 'Earth':
        
        timemsg = "Choose image time method"
        timechoices = ["User Entry","Yudish Imagetimeheader"]
        reply = easygui.buttonbox(timemsg, choices=timechoices)
        
        if reply == "Yudish Imagetimeheader": [idls,ctime,uncertainty_range_exists] = image_time_yudish(comdenom,fitsinfile,idlsav)
        if reply == "User Entry": [ctime,uncertainty_range_exists] = image_time_user()
                                          
    elif obsloc == 'Stereo_A' or 'Stereo_B':
        
        [ctime,uncertainty_range_exists] = image_time_stereo(filebase)

    elif obsloc == 'Soho': sys.exit("SOHO Image times need adding")
        
    #find relevant observer parameters of comet at observer time
    comcel = np.where(abs(obsveceq[:,0] - ctime.jd) < 1e-4)[0][0]
    
    #find closest matching cell in dt = 10min data
    comcel10 = np.where(abs(comveceq10[:,0] - ctime.jd) < 3.5e-3)[0][0]
    c10t = astropy.time.Time(comveceq10[comcel10,0], format = 'jd')
    
    #find distance to actual time
    dt = ctime - c10t
    dtmin = int(round(dt.sec/60))
    
    write_properties(d,border, pixwidth, fnt, featur_fill, obsloc, filebase,
                     comcel, ctime, comobs)
                     
#%%**********************************************************
#FIFTH CELL - Plot Comet Traj, Comet-Sun Vector + Uncertainty
#************************************************************

    #ALSO NEED TO ADD SYNCHRONE TESTING AS IN LINUX

    #find ra and dec of comet
    LT_cor = int(np.round(np.linalg.norm(comveceq[comcel,6:9] - 
    obsveceq[comcel,6:9])*8.316746397269274))
    com_ra_dec = pos2radec(comveceq[comcel - LT_cor,6:9] - obsveceq[comcel,6:9])
    
    #check if comet is within image
    com_box_path = np.array(
                    	[[ra_m[0,0],		dec[0,0]		],
                   		 [ra_m[ya-1,0],	 	dec[ya-1,0]		],
                  		 [ra_m[ya-1,xa-1],  dec[ya-1,xa-1]	],
                   		 [ra_m[0,xa-1],	 	dec[0,xa-1]		]])              		 
    com_Path = mplPath.Path(com_box_path)    
    
    com_in_image = com_Path.contains_point((com_ra_dec[0],com_ra_dec[1]))
    if com_in_image == True:
     
        #setting RGBA colours of lines
        trajfill = (0,255,255,255)
        trajucfill = (0,100,255,255)
        comsunfill = (0,255,0,255)
    
        #account for sometimes negative axisdata
        if (axisdata[6] < 0):
            ra_img_lower = axisdata[6] + 360
        else: ra_img_lower = axisdata[6]
        if (axisdata[7] < 0):
            ra_img_higher = axisdata[7] + 360
        else: ra_img_higher = axisdata[7]
        
        [ltcomcel, vtraj, vtrajcel] = plot_orbit(comobs,comveceq,obsveceq,axisdata,d,comcel,trajfill,ra_img_lower,
            ra_img_higher,border,pixwidth,rafmin,scale,pixheight,decmin,fnt,featur_fill)        

        plot_orbit_points(d,vtraj,largefnt,featur_fill,pixheight,pixwidth,border)
        
        plot_sunearth_vec(d,comveceq,obsveceq,ltcomcel,axisdata,ra_img_higher,ra_img_lower,
                        border,pixwidth,pixheight,rafmin,decmin,scale,comsunfill,featur_fill,fnt)
        
        if obsloc == 'Earth' and uncertainty_range_exists == True:
            plot_compos_unc(d,idls,vtraj,vtrajcel,border,pixwidth,fnt,trajucfill,featur_fill)
    else:
        trajfill = None;trajucfill = None;comsunfill = None
        
#%%***********************************************************
#SIXTH CELL - Save image and parameters or load existing image
#*************************************************************

    rapixl = ra2xpix(com_ra_dec[0],border,pixwidth,rafmin,scale)
    decpixl = dec2ypix(com_ra_dec[1],border,pixheight,decmin,scale)

    comimg.show() #shows for refrence
    comimg.save(imgsav,'png')
    
    with open(picklesavefile, 'w') as f:
        pickle.dump([comcel, comcel10, ramax, decmax, ramin, decmin, border,
                     pixheight, pixwidth, scale, ctime, dtmin, ra, dec, colr,
                     colb, colg, trajfill, trajucfill, comsunfill, backgr_fill,
                     imgwidth, imgheight, rafmin, rafmax, rapixl, decpixl,
                     com_ra_dec, com_Path, com_in_image, featur_fill], f)
                    
else:
    print "Loading parameter data"
    
    with open(picklesavefile) as f:
        parameters = pickle.load(f)
        comcel = parameters[0]; comcel10 = parameters[1]; ramax = parameters[2]
        decmax = parameters[3]; ramin = parameters[4]; decmin = parameters[5]
        border = parameters[6]; pixheight = parameters[7]; pixwidth = parameters[8]
        scale = parameters[9]; ctime = parameters[10]; dtmin = parameters[11]
        trajfill = parameters[17]; trajucfill = parameters[18]; comsunfill = parameters[19]
        backgr_fill = parameters[20]; imgwidth = parameters[21]; imgheight = parameters[22]
        rafmin = parameters[23]; rafmax = parameters[24]; rapixl = parameters[25]
        decpixl = parameters[26]; com_ra_dec = parameters[27]; com_Path = parameters[28]
        com_in_image = parameters[29]; featur_fill = parameters[30]
        
    comimg = Image.open(imgsav)
    comimg.show()
   
#%%**********************************************
#SEVENTH CELL - Preparing to simulate dust motion
#************************************************

test_mode = True
while test_mode == True:
    
    comimg = Image.open(imgsav)
    d = ImageDraw.Draw(comimg)
    
    simsavefile = os.path.join(pysav, filebase + '_simsetup')
    [betau, betal, bno, simtu, simtl, tno, tspace, drawopts, sav_bool, test_mode] \
    = simulation_setup(simsavefile)
    
    #finds maximum possible simulateable time and ensure simtu is bounded by this
    simtmax = np.round(ctime.jd - comveceq10[0,0])
    if (simtu > simtmax):
        simtu = simtmax
        print "Simulation time exceeds orbit data range, setting to max"
    else:
        simtu = simtu
        
    #option for log distributed or linear distributed
    if (tspace == 'Logarithmic'):
        tvals = np.logspace(np.log10(simtl), np.log10(simtu), num=tno)
    elif (tspace == 'Linear'):
        tvals = np.linspace(simtl, simtu, tno)
    
    #create log distributed beta values, accounting that log(0) is impossible
    if (betal == 0):
        bvals = np.logspace(np.log10(0.001), np.log10(betau), num=(bno-1))
        bvals = np.concatenate((np.array([0]), bvals))
        print "Beta = 0, logarithmic spacing from Beta = 0.001"
    else:
        bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))

    #prep for simulation loop
    simres = np.empty((tno,bno,16),dtype = float)
    efinp = obsveceq[comcel,6:9]
    bmax = np.empty((tno),dtype = int)
    tidx = 0
    
#%%*********************************
#EIGHT CELL - Simulating Dust Motion
#***********************************
    
    #this is the loop that does the business
    while (tidx < tno):
        bidx = 0
        point_in_image = 1
        while ( bidx < bno and point_in_image == 1):
            simt10min = int(round(144*tvals[tidx]))
            pstart = comveceq10[comcel10-simt10min,6:12]
            sim = part_sim(bvals[bidx],simt10min,30,3,pstart,efinp,dtmin)
            simres[tidx,bidx,0] = float(simt10min)/144
            simres[tidx,bidx,1] = bvals[bidx]
            simres[tidx,bidx,2] = sim[0] #length of simulation in minutes
            simres[tidx,bidx,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
            simres[tidx,bidx,4:10] = sim[1] #finishing pos/vel
            simres[tidx,bidx,10:12] = pos2radec(sim[1][0:3] - 
            obsveceq[int(simres[tidx,bidx,3]),6:9])
            point_in_image = int(com_Path.contains_point(
            (simres[tidx,bidx,10],simres[tidx,bidx,11])))
            simres[tidx,bidx,12] = ra2xpix(rasim,border,pixwidth,rafmin,scale)
            simres[tidx,bidx,13] = dec2ypix(desim,border,pixheight,decmin,scale)                 
            simres[tidx,bidx,14] = 1
            bidx += 1
        tidx += 1
        print float(tidx)*100/tno
        
    tmax = tno - 1 - np.searchsorted(bmax[::-1], np.arange(bno))

    if (sav_bool == True):
        simressavefile = os.path.join(pysav, 'simres')
        simressavefile = os.path.join(simressavefile, filebase + '_' + str(betal) + '_'
                         + str(betau) + '_' + str(bno)+ '_' + str(simtu) + '_'
                         + str(simtl) + '_' + str(tno))
        simressavefile = string.replace(simressavefile,'.','\'')
        np.save(simressavefile, simres)
        with open(simressavefile + '_parameters' , 'w') as f:
            pickle.dump([tmax, bmax, tno, bno, tspace], f)
    
#%%***************************
#NINTH CELL - Plot dust motion
#*****************************

    dynfill = (255,0,0,255) #syndynes
    chrfill = (255,192,0,255) #synchrones
    drfill = (255,0,255,255) #data points and data region
    
    fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
    fnt = ImageFont.truetype(fontloc, 20)  
    bt_anno_idx = 0
        
    if "Syndynes" in drawopts: draw_synchrones(dynfill,d,simres,bno,rapixl,decpixl,tmax)
    elif "Synchrones" in drawopts: draw_sydynes(chrfill,d,simres,tno,rapixl,decpixl,bmax)
    elif "Data Points" in drawopts: draw_datap(drfill,d,simres,tno,bmax)
    elif "Data Region Enclosed" in drawopts:  draw_data_reg(drfill,featur_fill,d,fnt,simres,tno,bmax,border,pixwidth)                
    bt_anno_idx = annotate_plotting(d,drawopts,border,pixwidth,fnt,featur_fill,dynfill,chrfill,drfill)
    
    write_bt_ranges(d,border, pixwidth, fnt, featur_fill,
                    betau, betal, simtu, simtl, bt_anno_idx)
    
    if (drawopts != "No Image"):
        comimg.show()    
        cimgsav = os.path.join(imagedir, 'temp.png')    
        comimg.save(cimgsav,'png')