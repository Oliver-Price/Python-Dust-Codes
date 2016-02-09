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
import sys
from astropy.io import fits
from astropy import wcs
import astropy.time
import datetime
from PIL import Image, ImageDraw, ImageFont
from orbitdata_loading_functions import orb_vector, orb_obs
from plot_functions import ra2xpix, dec2ypix, setaxisup
from conversion_routines import pos2radec
from Overplot_Simulation_Setup import simulation_setup
from particle_sim import part_sim
import idlsave
import pickle

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
    obsloc = cdata[34][19:24]
    horiztag = cdata[40][10:]

#import the orbit data
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,eq')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'eq,d10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

#choosing fits file to display and getting pathnames
fitsin = easygui.fileopenbox(default = os.path.join(imagedir, r"*.fits"))
fitsinfile = os.path.basename(fitsin)
filebase = fitsinfile[:string.find(fitsinfile,'.')]

inv = False

#check if image exists already
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
    
    #choose img type
    colormsg = "Please select image type:"
    colorchoices = ["Colour","Black & White","Quit"]
    reply = easygui.buttonbox(colormsg, choices=colorchoices)
        
    if reply == "Colour":

        #we have to use a temporary 1d fits file to get WCS data
        #this is due to astropy being unable to handle RGB fits images well
        fitstemp = os.path.join(imagedir, 'temporary_' + fitsinfile)
        hdulist = fits.open(fitsin)
        
        #opening color data
        colours = (hdulist[0].data)
        colr = colours[0,:,:]
        colg = colours[1,:,:]
        colb = colours[2,:,:]
        
        #making a 1d image for doing ra/dec coordinates
        oldtemp = os.path.isfile(fitstemp) #if one doesn't already exist
        if oldtemp == False: 
            hdulist[0].data = hdulist[0].data[0] #makes image from red plane
            hdulist[0].header['naxis'] = 1       #edits header so ds9 will work
            hdulist.writeto(fitstemp)
            
        hdulist.close()
        fitscoords = fitstemp #directs program to look at temp image for coords
     
    elif reply == "Black & White":
        
        #simple case, we can use base image for coords
        hdulist = fits.open(fitsin)
        colours = (hdulist[0].data)
        fitscoords = fitsin
        colr = colours #black and white colour scheme
        colg = colours
        colb = colours
    
    #otherwise, stop the program
    elif reply == "Quit":
        sys.exit()
        
    if inv == False:   
        colcr = colr; colcg = colg; colcb = colb
        backgr_fill = (0,0,0,255)
        featur_fill = (255,255,255,255)
    elif inv == True:
        colcr = 255 - colr; colcg = 255 - colg; colcb = 255 - colb
        backgr_fill = (255,255,255,255)
        featur_fill = (0,0,0,255)
        
#%%**********************
#SECOND CELL - Plot Image
#************************

    #get RA/DEC data    
    onedimg = fits.open(fitscoords)
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
    radecs = w.wcs_pix2world(coords, 1)
    ra = np.reshape(radecs[:,0], (ya,xa))
    dec = np.reshape(radecs[:,1], (ya,xa))
    
    #find minimum/maximum values
    ramin = np.amin(ra)
    ramax = np.amax(ra)
    decmin = np.amin(dec)
    decmax = np.amax(dec)
    
    #make a canvas with a fixed pixel height and border
    pixheight = 800
    pixwidth = int(pixheight*(ramax - ramin)/(decmax - decmin))
    border = 100
    scale = pixheight/(decmax - decmin)
    imgwidth = pixwidth+int(4*border)
    imgheight = pixheight+int(3*border)
    comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,backgr_fill)
    d = ImageDraw.Draw(comimg)
    
    #plots image on canvas
    for x in xrange(0, ya-2):
        for y in xrange(0, xa-2):
            ra1 = ra2xpix(ra[x,y],border,pixwidth,ramin,scale)
            ra2 = ra2xpix(ra[x+1,y],border,pixwidth,ramin,scale)
            ra3 = ra2xpix(ra[x+1,y+1],border,pixwidth,ramin,scale)
            ra4 = ra2xpix(ra[x,y+1],border,pixwidth,ramin,scale)
            dec1 = dec2ypix(dec[x,y],border,pixheight,decmin,scale)
            dec2 = dec2ypix(dec[x+1,y],border,pixheight,decmin,scale)
            dec3 = dec2ypix(dec[x+1,y+1],border,pixheight,decmin,scale)
            dec4 = dec2ypix(dec[x,y+1],border,pixheight,decmin,scale)
            a = d.polygon([(ra1,dec1),(ra2,dec2),(ra3,dec3),(ra4,dec4)] ,\
            fill=(colcr[x,y],colcg[x,y],colcb[x,y],255))
            
    #%%********************
    #THIRD CELL - Draw Axis
    #**********************
    
    #draws a border       
    a = d.polygon([(border,border),(border*2+pixwidth,border), \
    (border*2+pixwidth,border*2+pixheight),(border,border*2+pixheight)], \
    outline = featur_fill)
    
    #most of the dirty stuff is bunged into this function
    axisdata = setaxisup(ramax,ramin,decmax,decmin,border,pixheight,pixwidth,scale)
    
    #axisdata 0/1 - ra major divisions --- values and pixel locations
    #axisdata 2 - ra minor divisions --- pixel locations
    #axisdata 3/4 - dec major divisions --- values and pixel locations
    #axisdata 5 - dec minor divisions --- pixel locations
    #axisdata 6/7 min/max RA in border
    #axisdata 8/9 min/max DEC in border
    
    majt = 20  #major tick length
    mint = 10  #minor tick length
    rdaxis = pixheight + border*2
    
    fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
    fnt = ImageFont.truetype(fontloc, 20)
    
    for div in xrange(0, (np.size(axisdata[1]))): #RA axis major ticks
        b = d.line([(axisdata[1][div],rdaxis-majt),(axisdata[1][div],rdaxis)],\
        fill = featur_fill)
        tick = str(axisdata[0][div])
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
    tfnt = ImageFont.truetype(fontloc, 30)
    d.text((1.5*border + pixwidth*0.5 - len(plttitle)*5 - 60,.35*border), \
    plttitle, font=tfnt, fill= featur_fill)
    
    pix = comimg.load()
    
    #%%**************************
    #FOURTH CELL - Get Imgtimehdr
    #****************************
    
    #locate imgtimehdr save data from IDL
    idlsavpath = os.path.join(idlsav,'Imagetimeheaders_savefile')
    idlsavpath = os.path.join(idlsavpath, comdenom)
    idlsavpath = os.path.join(idlsavpath, fitsinfile[len(comdenom)+1:-5])
    idlsavpath = idlsavpath + '_timeinfo.sav'
    idls = idlsave.read(idlsavpath)
    
    #get time of comet in image and make an astropy time reference
    chour = int(idls.optocentre_time_str[0][0:2])
    cmin = int(idls.optocentre_time_str[0][3:5])
    cday = int(idls.optocentre_date[0][0:2])
    cmonth = int(idls.optocentre_date[0][3:5])
    cyear = int(idls.optocentre_date[0][6:11])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                chour , cmin, 0))
    
    #find relevant observer parameters of comet at observer time
    comcel = np.where(abs(obsveceq[:,0] - ctime.jd) < 1e-4)[0][0]
    
    #find closest matching cell in dt = 10min data
    comcel10 = np.where(abs(comveceq10[:,0] - ctime.jd) < 3.5e-3)[0][0]
    c10t = astropy.time.Time(comveceq10[comcel10,0], format = 'jd')
    
    #find distance to actual time
    dt = ctime - c10t
    dtmin = int(round(dt.sec/60))
    
    #display orbit plane angle
    cimgopa = comobs[comcel,7]
    d.text((2*border + pixwidth + 30,border), \
    "Plane Angle:" + '\n' + '  ' + "%.2f" % cimgopa ,
    font=fnt, fill= featur_fill)
    
    #display image author
    d.text((2*border + pixwidth + 30,border + 50), \
    "Image From: " + '\n' + '  ' + filebase.split('_')[-1] ,
    font=fnt, fill= featur_fill)
    
    #display date
    d.text((2*border + pixwidth + 30,border + 100), \
    "Image Date: " + '\n' + '  ' + ctime.isot[0:10] ,
    font=fnt, fill= featur_fill)

    #display time
    d.text((2*border + pixwidth + 30,border + 150), \
    "Image Time: " + '\n' + '  ' + ctime.isot[11:16] ,
    font=fnt, fill= featur_fill)    

    #%%**********************************************************
    #FIFTH CELL - Plot Comet Traj, Comet-Sun Vector + Uncertainty
    #************************************************************
    
    #setting RGBA colours of lines
    trajfill = (0,255,255,255)
    trajucfill = (0,100,255,255)
    comsunfill = (0,255,0,255)
    
    #find rough cell range of traj from observer data
    obsraloc = np.where((comobs[:,5] > axisdata[6]) \
    & (comobs[:,5] < axisdata[7]))[0]
    obsdecloc = np.where((comobs[:,6] > axisdata[8]) \
    & (comobs[:,6] < axisdata[9]))[0]
    trajrough = np.intersect1d(obsraloc,obsdecloc)
    
    #use this to calculate ra and dec of comet for a purposely oversized range
    vno = 0; vext = int(np.size(trajrough))
    vtraj = np.empty((np.size(trajrough)+2*vext-1,6),dtype = float)
    tcellmax = min(trajrough[-1] + vext, np.shape(comveceq)[0])
    for tcell in xrange(trajrough[0] - vext,tcellmax):
        vtemp = comveceq[tcell,6:9] - obsveceq[comcel,6:9]    
        ptemp = pos2radec(vtemp)
        vtraj[vno,0] = ptemp[0]
        vtraj[vno,1] = ptemp[1]
        vtraj[vno,4] = tcell
        vtraj[vno,5] = tcell + round(np.linalg.norm(vtemp)*8.316746397269274)
        vno +=1
        
    #use these ra and dec values to slim data down to within image borders
    trajrange = np.intersect1d(
                np.intersect1d( np.where(vtraj[:,0] < axisdata[7])[0],
                                np.where(vtraj[:,0] > axisdata[6])[0]),
                np.intersect1d( np.where(vtraj[:,1] < axisdata[9])[0],
                                np.where(vtraj[:,1] > axisdata[8])[0]))                  
    vtraj = vtraj[trajrange[0]:trajrange[-1],:]
    
    #find relevant cell in vtraj and comveceq accounting for LT
    vtrajcel = np.where(abs(vtraj[:,5] - comcel) < 1e-4)[0][0]
    ltcomcel = vtraj[vtrajcel,4]
    
    #convert to ra and dec, and plot
    vtraj[:,2] = ra2xpix(vtraj[:,0],border,pixwidth,ramin,scale)
    vtraj[:,3] = dec2ypix(vtraj[:,1],border,pixheight,decmin,scale)
    for ta in xrange(0, (np.shape(vtraj)[0]-1)):
        b = d.line([(vtraj[ta,2],vtraj[ta,3]),(vtraj[ta+1,2],vtraj[ta+1,3])],\
        fill = trajfill)
    
    #creates an array of ra/dec values along sun-comet line
    cmsam = 10001
    offsets = np.ones(cmsam) + np.linspace(-0.3,0.7,cmsam)
    comsun = np.empty((cmsam,4),dtype = float)
    for cs in xrange(0,cmsam):
        ctemp = pos2radec(offsets[cs]*comveceq[ltcomcel,6:9] - 
                            obsveceq[ltcomcel,6:9])
        comsun[cs,0] = ctemp[0]
        comsun[cs,1] = ctemp[1]
    
    #slims this array down to within image limits
    csvrange = np.intersect1d(
               np.intersect1d( np.where(comsun[:,0] < axisdata[7])[0],
                                np.where(comsun[:,0] > axisdata[6])[0]),
               np.intersect1d( np.where(comsun[:,1] < axisdata[9])[0],
                                np.where(comsun[:,1] > axisdata[8])[0]))                  
    comsun = comsun[csvrange[0]:csvrange[-1],:]
    
    #convert to ra and dec, and plot
    comsun[:,2] = ra2xpix(comsun[:,0],border,pixwidth,ramin,scale)
    comsun[:,3] = dec2ypix(comsun[:,1],border,pixheight,decmin,scale)
    for ca in xrange(0, (np.shape(comsun)[0]-1)):
        b = d.line([(comsun[ca,2],comsun[ca,3]),
                    (comsun[ca+1,2],comsun[ca+1,3])],
                    fill = comsunfill)
        
    #draw range of uncertainty in comet position
    sig_t = int(float(idls.sig_t))
    ura = vtraj [ (vtrajcel - sig_t) : (vtrajcel + sig_t) , 2 ]
    udec = vtraj [ (vtrajcel - sig_t) : (vtrajcel + sig_t) , 3 ]
    for ua in xrange(0, np.size(ura)-2):
        b = d.line([(ura[ua],udec[ua]),(ura[ua+1],udec[ua+1])],  
        fill = trajucfill)
        
    d.text((2*border + pixwidth + 30,border + 200), \
    "Orbital Path:",font=fnt, fill= featur_fill)
    d.line([(2*border + pixwidth + 30,border + 230),
            (2*border + pixwidth + 170,border + 230)], fill = trajfill)
            
    d.text((2*border + pixwidth + 30,border + 250), \
    "Orbital\nUncertainty:",font=fnt, fill= featur_fill)
    d.line([(2*border + pixwidth + 30,border + 294),
            (2*border + pixwidth + 170,border + 294)], fill = trajucfill)
            
    d.text((2*border + pixwidth + 30,border + 310), \
    "Comet-Sun\nVector:",font=fnt, fill= featur_fill)
    d.line([(2*border + pixwidth + 30,border + 355),
            (2*border + pixwidth + 170,border + 355)], fill = comsunfill) 
            
            
#%%***********************************************************
#SIXTH CELL - Save image and parameters or load existing image
#*************************************************************
    
    #store parameters from vtraj used later on as floats to save space
    rao = vtraj[vtrajcel,0]
    deo = vtraj[vtrajcel,1]
    rapixl = vtraj[vtrajcel,2]
    depixl = vtraj[vtrajcel,3]

    comimg.show() #shows for refrence
    comimg.save(imgsav,'png')
    
    with open(picklesavefile, 'w') as f:
        pickle.dump([comcel, comcel10, rao, deo, rapixl, depixl, ramax, decmax
                    ,ramin, decmin, border, pixheight, pixwidth, scale, ctime
                    ,dtmin, ra, dec, colr, colb, colg, trajfill, trajucfill,
                    comsunfill, backgr_fill, imgwidth, imgheight], f)
                    
else:
    print "Loading save image and parameter data"
    
    with open(picklesavefile) as f:
        parameters = pickle.load(f)
        comcel = parameters[0]
        comcel10 = parameters[1]
        rao = parameters[2]
        deo = parameters[3]
        rapixl = parameters[4]
        depixl = parameters[5]
        ramax = parameters[6]
        decmax = parameters[7]
        ramin = parameters[8]
        decmin = parameters[9]
        border = parameters[10]
        pixheight = parameters[11]
        pixwidth = parameters[12]
        scale = parameters[13]
        ctime = parameters[14]
        dtmin = parameters[15]
        trajfill = parameters[21]
        trajucfill = parameters[22]
        comsunfill = parameters[23]
        backgr_fill = parameters[24]
        imgwidth = parameters[25]
        imghegiht = parameters[26]
        
    comimg = Image.open(imgsav)
    d = ImageDraw.Draw(comimg)
    comimg.show()
   
#%%**********************************************
#SEVENTH CELL - Preparing to simulate dust motion
#************************************************

simsavefile = os.path.join(pysav, filebase + '_simsetup')
[betau, betal, bno, simtu, simtl, tno, tspace, threshold, drawopts] = \
simulation_setup(simsavefile)

#finds maximum possible simulateable time and ensure simtu is bounded by this
simtmax = np.round(ctime.jd - comveceq10[0,0])
if (simtu > simtmax):
    simtu = simtmax 
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
else:
    bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))

#this is for identifying pixels where lines are present on grid
linevalfalses = np.ones((256,256,256), dtype = int)
linevalfalses[trajfill[0],trajfill[1],trajfill[2]] = 0
linevalfalses[trajucfill[0],trajucfill[1],trajucfill[2]] = 0
linevalfalses[comsunfill[0],comsunfill[1],comsunfill[2]] = 0
linevalfalses[backgr_fill[0],backgr_fill[1],backgr_fill[2]] = 0

#prep for simulation loop
simres = np.empty((tno,bno,16),dtype = float)
efinp = obsveceq[comcel,6:9]
bmax = np.empty((tno),dtype = int)
tidx = 0
pix = comimg.load()

#%%*********************************
#EIGHT CELL - Simulating Dust Motion
#***********************************

while (tidx < tno):
    bidx = 0
    rasim = rao
    desim = deo
    while ( bidx < bno and rasim <= ramax and rasim >= ramin  and
            desim <= decmax and desim >= decmin):
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
        rasim = simres[tidx,bidx,10]
        desim = simres[tidx,bidx,11]
        simres[tidx,bidx,12] = ra2xpix(rasim,border,pixwidth,ramin,scale)
        simres[tidx,bidx,13] = dec2ypix(desim,border,pixheight,decmin,scale)
        pv = pix[ int( round( sorted([0, simres[tidx,bidx,12], pixwidth])[1] ) ),
        int( round(  sorted([0, simres[tidx,bidx,13], pixheight])[1] ) ) ][0:3]        
        simres[tidx,bidx,14] = linevalfalses[pv[0],pv[1],pv[2]]*np.average(pv)                     
        simres[tidx,bidx,15] = 1
        bidx += 1
    simres[tidx,bidx-1,15] = 0
    bmax[tidx] = bidx - 1
    tidx += 1
    print float(tidx)*100/tno
    
tmax = tno - 1 - np.searchsorted(bmax[::-1], np.arange(bno))

if inv == True:
        zero_locs = np.where(simres[:,:,14] <= threshold)
elif inv == False:
        zero_locs = np.where(simres[:,:,14] >= (255 - threshold))

zero_size = np.size(zero_locs[0])
for z in xrange(0, zero_size):
    simres[zero_locs[0][z],zero_locs[1][z],15] = 0

sav = True
if (sav == True):
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

if inv ==  False:
    sfill = (255,0,0,255)
elif inv ==  True:
    sfill = (255,0,0,255)
    
if "Syndynes" in drawopts: #DRAW SYNDYNES
    for ba in xrange(0, bno):
        b = d.line([(rapixl,depixl), \
        (simres[0,ba,12],simres[0,ba,13])],\
        fill = sfill)
        for ta in xrange(0, tmax[ba]):
            b = d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
            (simres[ta+1,ba,12],simres[ta+1,ba,13])],\
            fill = sfill)

if "Synchrones" in drawopts: #DRAW SYNCHRONES
    for ta in xrange(0, tno):
        b = d.line([(rapixl,depixl), \
        (simres[ta,0,12],simres[ta,0,13])],\
        fill = sfill)
        for ba in xrange(0, bmax[ta]):
            b = d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
            (simres[ta,ba+1,12],simres[ta,ba+1,13])],\
            fill = sfill)

if "Data Points" in drawopts: #DRAW DATAPOINTS
    for ta in xrange(0,tno):
        for ba in xrange(0, bmax[ta]):
            xsiz = 2
            b = d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] - xsiz ) ,
                      ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] + xsiz ) ] ,
                      fill = (255,0,255,255) )  
            b = d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] + xsiz ) ,
                      ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] - xsiz ) ] ,
                      fill = (255,0,255,255) )

if drawopts == "Data Region Enclosed":
    b = d.line( [ ( simres[0,0,12]  , simres[0,0,13] ) ,
                  ( simres[tno-1,0,12] , simres[tno-1,0,13] ) ] ,
                    fill = (255,0,255,255) )
    b = d.line( [ ( simres[0,0,12]  , simres[0,0,13] ) ,
                  ( simres[0,bmax[0],12] , simres[0,bmax[0],13] ) ] ,
                    fill = (255,0,255,255) )
    b = d.line( [ ( simres[tno-1,0,12]  , simres[tno-1,0,13] ) ,
                  ( simres[tno-1,bmax[tno-1],12] ,
                   simres[tno-1,bmax[tno-1],13] ) ] ,    
                    fill = (255,0,255,255) )           
    for ta in xrange(0,tno - 1):
        b = d.line( [ ( simres[ta,bmax[ta],12]  , simres[ta,bmax[ta],13] ) ,
                      ( simres[ta+1,bmax[ta+1],12] , simres[ta+1,bmax[ta+1],13] ) ]
                      , fill = (255,0,255,255) )

if (drawopts != "No Image"):
    comimg.show()    
    cimgsav = os.path.join(imagedir, 'temp.png')    
    comimg.save(cimgsav,'png')