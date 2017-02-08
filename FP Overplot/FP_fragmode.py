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
import webbrowser
import matplotlib.pyplot as plt

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from orbitdata_loading_functions import orb_vector, orb_obs
from FP_plot_functions import ra2xpix, dec2ypix, setaxisup, plotpixel, getdustphase
from FP_plot_functions import draw_syndynes, draw_synchrones, draw_datap
from FP_plot_functions import draw_data_reg, annotate_plotting
from FP_plot_functions import draw_phase_points, annotate_dustphase
from FP_diagnostics import plot_orbit, plot_orbit_points, plot_sunearth_vec
from FP_diagnostics import write_bt_ranges, write_properties, plot_compos_unc
from imagetime_methods import image_time_yudish, image_time_user, image_time_stereo, image_time_filename
from conversion_routines import pos2radec, fixwraps, find_largest_nonzero_block 
from io_methods import get_obs_loc, correct_for_imagetype, get_hih_low
from io_methods import get_stereo_instrument, get_soho_instrument
from simulation_setup import simulation_setup
from particle_sim import part_sim, frag_sim

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

#choosing fits file to display and getting pathnames
fitsin = easygui.fileopenbox(default = os.path.join(imagedir,'*'))
if fitsin == '.': sys.exit("No file selected")
fitsinfile = os.path.basename(fitsin)
filebase = fitsinfile[:string.find(fitsinfile,'.')]

inv = False
pngdir = os.path.join(imagedir, 'cometplots')
if not os.path.exists(pngdir): os.makedirs(pngdir)
#check if image exists already/invert it
if inv == False:   
    imgsav = os.path.join(pngdir, filebase + '.png')
    imgexists = os.path.exists(imgsav)
elif inv == True:
    imgsav = os.path.join(pngdir, filebase + '_inverted.png')
    imgexists = os.path.exists(imgsav)  

#parameter savefile locations
picklesavefile = os.path.join(pysav, 'imgsavs')
if not os.path.exists(picklesavefile): os.makedirs(picklesavefile)
picklesavefile = os.path.join(picklesavefile, obsloc)
if not os.path.exists(picklesavefile): os.makedirs(picklesavefile)
if "Stereo" in obsloc: picklesavefile = os.path.join(picklesavefile, inst)
if "Soho" in obsloc: picklesavefile = os.path.join(picklesavefile, inst)
if not os.path.exists(picklesavefile): os.makedirs(picklesavefile)
picklesavefile = os.path.join(picklesavefile, filebase + '_plot_param.pickle')
 
picklexists = os.path.exists(picklesavefile)

forceredraw = False
if (imgexists == False) or (picklexists == False) or (forceredraw == True):
    print "preparing image"
    
    #ensures image inputted correctly depending on size of data cube
    [colr, colg, colb, fitscoords] = correct_for_imagetype(imagedir, fitsin, fitsinfile)
    
    #correcting for wrong data type, account for inversion FOR NON FITS ONLY
    if inv == False:
        backgr_fill = (0,0,0,255)
        featur_fill = (255,255,255,255)
        trajfill = (0,255,255,255)
        trajucfill = (0,0,255,255)
    elif inv == True:
        backgr_fill = (255,255,255,255)
        featur_fill = (0,0,0,255)
        trajfill = (0,0,255,255)
        trajucfill = (0,255,255,255)
        
#%%**********************
#SECOND CELL - Plot Image
#************************

    #get RA/DEC data    
    onedimg = fits.open(fitscoords)

    if 'Stereo' in obsloc:
        if comdenom == 'c2011l4':
            w = wcs.WCS(onedimg[0].header)
        elif comdenom == 'c2006p1':
            w = wcs.WCS(onedimg[0].header, key = 'A')
        if 'diff' in inst or 'MGN' in inst: plotmethodlog = False
        else: plotmethodlog = True
    elif 'Soho' in obsloc:
        w = wcs.WCS(onedimg[0].header)
        if 'diff' in inst or 'MGN' in inst: plotmethodlog = False
        else: plotmethodlog = True
    elif 'Earth' in obsloc:
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
    comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,backgr_fill)
    d = ImageDraw.Draw(comimg)
    
    if ("Stereo" in obsloc) and ('2' in inst):
        maskedimg = fits.open('C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\hi2A_mask.fts')
        imagemask  = maskedimg[0].data[::2,::2][:-1,:-1]
    else: imagemask = np.ones_like(colr)[:-1,:-1]
    
    [hih,low] = get_hih_low(comdenom,obsloc,inst)
    
    if plotmethodlog == True:
        
        grad = (255 / (np.log10(hih) - np.log10(low)))      
        colcr = (np.clip(np.round((np.log10(np.clip(colr,1e-20,999999999))- np.log10(low))*grad),0,255)).astype(int)
        colcg = (np.clip(np.round((np.log10(np.clip(colg,1e-20,999999999))- np.log10(low))*grad),0,255)).astype(int)
        colcb = (np.clip(np.round((np.log10(np.clip(colb,1e-20,999999999))- np.log10(low))*grad),0,255)).astype(int)
        non_zeros_0 = np.where(imagemask == 1)[0]
        non_zeros_1 = np.where(imagemask == 1)[1]
        no_points = non_zeros_0.size
        
        if inv == True:
            colcr = 255-colcb; colcb = 255-colcb; colcg = 255-colcg
        
        for xp in xrange(0,no_points):
            x = non_zeros_0[xp]
            y = non_zeros_1[xp]
            plotpixel(d,x,y,ra_m,dec,border,pixwidth,pixheight,decmin,rafmin,
                          scale,colcr[x,y],colcg[x,y],colcb[x,y])
                
    elif plotmethodlog == False:

        colcr = np.clip(255.0/(hih-low)*(colr-low),0,255).astype(int)
        colcg = np.clip(255.0/(hih-low)*(colg-low),0,255).astype(int)
        colcb = np.clip(255.0/(hih-low)*(colb-low),0,255).astype(int)         
            
        non_zeros_0 = np.where(imagemask == 1)[0]
        non_zeros_1 = np.where(imagemask == 1)[1]
        no_points = non_zeros_0.size
        
        if inv == True:
            colcr = 255-colcb; colcb = 255-colcb; colcg = 255-colcg
        
        for xp in xrange(0,no_points):
            x = non_zeros_0[xp]
            y = non_zeros_1[xp]
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
    fnt = ImageFont.truetype(fontloc, 25)
    smallfnt = ImageFont.truetype(fontloc, 20)
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
        d.text((border - 8 - 15*len(tick),axisdata[4][div] - 10 ), \
        tick, font=fnt, fill=featur_fill)
        
    for div in xrange(0, (np.size(axisdata[5]))): #DEC axis minor ticks
        b = d.line([(border+mint,axisdata[5][div]),(border,axisdata[5][div])],\
        fill= featur_fill)
    
    #axis labels
    d.text((1.5*border + pixwidth*0.5 - 145,pixheight + 2.5*border), \
    "Right Ascension (Degrees)", font=fnt, fill= featur_fill)
    d.text((0.25*border - 10,0.75*border - 20), \
    "Declination (Degrees)", font=fnt, fill= featur_fill)
    '''
    #plot title
    plttitle = (comdenom.upper() + ' ' + comname)
    d.text((1.5*border + pixwidth*0.5 - len(plttitle)*5 - 60,.35*border), \
    plttitle, font=largefnt, fill= featur_fill)
    '''
    pix = comimg.load()
    
#%%**************************
#FOURTH CELL - Get Imgtimehdr
#****************************
    
    #some amateur images lack times, need to sometimes use yudish imgtimehdr
    if obsloc == 'Earth':
        
        timemsg = "Choose image time method"
        timechoices = ["Filename","User Entry","Yudish Imagetimeheader"]
        reply = easygui.buttonbox(timemsg, choices=timechoices)
        
        if reply == "Filename": [ctime,uncertainty_range_exists] = image_time_filename(filebase)
        if reply == "Yudish Imagetimeheader": [idls,ctime,uncertainty_range_exists] = image_time_yudish(comdenom,fitsinfile,idlsav)
        if reply == "User Entry": [ctime,uncertainty_range_exists] = image_time_user()
                                          
    elif 'S' in obsloc:
        
        [ctime,uncertainty_range_exists] = image_time_stereo(filebase)

    #find relevant observer parameters of comet at observer time
    comcel = np.where(abs(obsveceq[:,0] - ctime.jd)==abs(obsveceq[:,0] - ctime.jd).min())[0][0]
    
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
    
    [ltcomcel, vtraj, vtrajcel] = plot_orbit(comobs,comveceq,obsveceq,axisdata,d,comcel,trajfill,ra_img_lower,
    ra_img_higher,border,pixwidth,rafmin,scale,pixheight,decmin,fnt,featur_fill, com_in_image,fixwrapsbool)      
    
    if vtraj is not None:
        plot_orbit_points(d,vtraj,smallfnt,featur_fill,pixheight,pixwidth,border)    
    
    if com_in_image == True:
     
        #setting RGBA colours of lines
        comsunfill = (0,255,0,255)
        
        plot_sunearth_vec(d,comveceq,obsveceq,ltcomcel,axisdata,ra_img_higher,ra_img_lower,
                        border,pixwidth,pixheight,rafmin,decmin,scale,comsunfill,featur_fill,fnt,fixwrapsbool)
        
        if obsloc == 'Earth' and uncertainty_range_exists == True:
            plot_compos_unc(d,idls,vtraj,vtrajcel,border,pixwidth,fnt,trajucfill,featur_fill)
    
    else:
        trajucfill = None;comsunfill = None   
    
    #draw the sun if it's in a good place
    sun_pos = pos2radec(-obsveceq[comcel,6:9],fixwrapsbool)
    if sun_pos[0] <= rafmax and sun_pos[0] >= rafmin:
        if sun_pos[1] <= decmax and sun_pos[1] >= decmin:
            sun_x = ra2xpix(sun_pos[0],border,pixwidth,rafmin,scale)
            sun_y = dec2ypix(sun_pos[1],border,pixheight,decmin,scale)
            d.ellipse((sun_x-2, sun_y-2, sun_x+2, sun_y+2), fill = (255,255,0,255))
#%%***********************************************************
#SIXTH CELL - Save image and parameters or load existing image
#*************************************************************

    rapixl = ra2xpix(com_ra_dec[0],border,pixwidth,rafmin,scale)
    decpixl = dec2ypix(com_ra_dec[1],border,pixheight,decmin,scale)

    comimg.save(imgsav,'png')
    webbrowser.open(imgsav)
    
    with open(picklesavefile, 'w') as f:
        pickle.dump([comcel, comcel10, ramax, decmax, ramin, decmin, border,
                     pixheight, pixwidth, scale, ctime, dtmin, ra, dec, colr,
                     colb, colg, trajfill, trajucfill, comsunfill, backgr_fill,
                     imgwidth, imgheight, rafmin, rafmax, rapixl, decpixl,
                     com_ra_dec, com_Path, com_in_image, featur_fill, fitscoords], f)
                    
else:
    print "Loading parameter data"
    
    with open(picklesavefile) as f:
        parameters = pickle.load(f)
        comcel = parameters[0]; comcel10 = parameters[1]; ramax = parameters[2]
        decmax = parameters[3]; ramin = parameters[4]; decmin = parameters[5]
        border = parameters[6]; pixheight = parameters[7]; pixwidth = parameters[8]
        scale = parameters[9]; ctime = parameters[10]; dtmin = parameters[11]; ra = parameters[12]
        trajfill = parameters[17]; trajucfill = parameters[18]; comsunfill = parameters[19]
        backgr_fill = parameters[20]; imgwidth = parameters[21]; imgheight = parameters[22]
        rafmin = parameters[23]; rafmax = parameters[24]; rapixl = parameters[25]
        decpixl = parameters[26]; com_ra_dec = parameters[27]; com_Path = parameters[28]
        com_in_image = parameters[29]; featur_fill = parameters[30]; fitscoords = parameters[31]
        
    comimg = Image.open(imgsav)
    webbrowser.open(imgsav)
   
#%%**********************************************
#SEVENTH CELL - Preparing to simulate dust motion
#************************************************
comimg = Image.open(imgsav)
d = ImageDraw.Draw(comimg)

[ra_m, rafmin, rafmax, fixwrapsbool] = fixwraps(ra, ramax, ramin)    

#finds maximum possible simulateable time and ensure simtu is bounded by this
simt = 6
betastart = 0.5
n = 4
m = 25
bvals = np.logspace(np.log10(0.3), np.log10(3.0), num=(m))

ftvals = np.linspace(1,3,n)
bfvals = np.logspace(np.log10(1), np.log10(5), num=(n))

efinp = obsveceq[comcel,6:9]

#prep for simulation loop
simres = np.zeros((n,n,15),dtype = float)
simt10min = int(round(144*simt))
#%%*********************************
#EIGHT CELL - Simulating Dust Motion
#***********************************
for ft in xrange(n):
    for bf in xrange(n):
        pstart = comveceq10[comcel10-simt10min,6:12]
        sim = frag_sim(betastart,bfvals[bf],simt10min,ftvals[ft],30,3,pstart,efinp,dtmin)
        simres[ft,bf,0] = ftvals[ft]
        simres[ft,bf,1] = bfvals[bf]
        simres[ft,bf,2] = sim[0] #length of simulation in minutes
        simres[ft,bf,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
        simres[ft,bf,4:10] = sim[1] #finishing pos/vel
        simres[ft,bf,10:12] = pos2radec(sim[1][0:3] - 
        obsveceq[int(simres[ft,bf,3]),6:9],fixwrapsbool)
        point_in_image = int(com_Path.contains_point(
        (simres[ft,bf,10],simres[ft,bf,11])))
        simres[ft,bf,12] = ra2xpix(simres[ft,bf,10],border,pixwidth,rafmin,scale)
        simres[ft,bf,13] = dec2ypix(simres[ft,bf,11],border,pixheight,decmin,scale)                 
        simres[ft,bf,14] = point_in_image
    print float(ft)*5
simres[np.where(simres[:,:,14]==0)] = 0

sync_orig = np.zeros((m,15),dtype = float)
for b in xrange(m):
    pstart = comveceq10[comcel10-simt10min,6:12]
    sim = part_sim(bvals[b],simt10min,30,3,pstart,efinp,dtmin)
    sync_orig[b,0] = float(simt10min)/144
    sync_orig[b,1] = bvals[b]
    sync_orig[b,2] = sim[0] #length of simulation in minutes
    sync_orig[b,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
    sync_orig[b,4:10] = sim[1] #finishing pos/vel
    sync_orig[b,10:12] = pos2radec(sim[1][0:3] - 
    obsveceq[int(simres[ft,bf,3]),6:9],fixwrapsbool)
    point_in_image = int(com_Path.contains_point(
    (sync_orig[b,10],sync_orig[b,11])))
    sync_orig[b,12] = ra2xpix(sync_orig[b,10],border,pixwidth,rafmin,scale)
    sync_orig[b,13] = dec2ypix(sync_orig[b,11],border,pixheight,decmin,scale)                 
    sync_orig[b,14] = point_in_image

point_orig = np.zeros((15),dtype = float)
pstart = comveceq10[comcel10-simt10min,6:12]
sim = part_sim(betastart,simt10min,30,3,pstart,efinp,dtmin)
point_orig[0] = float(simt10min)/144
point_orig[1] = betastart
point_orig[2] = sim[0] #length of simulation in minutes
point_orig[3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
point_orig[4:10] = sim[1] #finishing pos/vel
point_orig[10:12] = pos2radec(sim[1][0:3] - 
obsveceq[int(simres[ft,bf,3]),6:9],fixwrapsbool)
point_in_image = int(com_Path.contains_point(
(point_orig[10],point_orig[11])))
point_orig[12] = ra2xpix(point_orig[10],border,pixwidth,rafmin,scale)
point_orig[13] = dec2ypix(point_orig[11],border,pixheight,decmin,scale)                 
point_orig[14] = point_in_image

#%%***************************
#NINTH CELL - Plot dust motion
#*****************************

xsiz = 5
d.line( [ ( point_orig[12] - xsiz , point_orig[13] - xsiz ) ,
          ( point_orig[12] + xsiz , point_orig[13] + xsiz ) ] ,
          fill = (255,150,50,255), width = 3 )  
d.line( [ ( point_orig[12] - xsiz , point_orig[13] + xsiz ) ,
          ( point_orig[12] + xsiz , point_orig[13] - xsiz ) ] ,
          fill = (255,150,50,255), width = 3 )

for ba in xrange(m-1):
    d.line([(sync_orig[ba,12],sync_orig[ba,13]), \
    (sync_orig[ba+1,12],sync_orig[ba+1,13])],\
    fill = (255,150,50,255),width = 2)
    
for ta in xrange(n):
    for ba in xrange(n-1):
        d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
        (simres[ta,ba+1,12],simres[ta,ba+1,13])],\
        fill = (int(255.*ta/(n-1)),255,0,255),width = 2)

fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
fnt = ImageFont.truetype(fontloc, 20)  
bt_anno_idx = 0

comimg.show()