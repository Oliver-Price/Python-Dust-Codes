#*************************************************************
#Program to visualise fits image of comet and overlay a
#finson-probstein diagram according to it's orbital parameters
#*************************************************************

import easygui
import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
import astropy.time
import pickle
import matplotlib.path as mplPath
from PIL import Image, ImageDraw, ImageFont
import webbrowser
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ollie\Documents\PhD\Python-Dust-Codes-Legacy\General-Use")
sys.path.append(r"C:\Users\Ollie\Documents\PhD\Python-Dust-Codes-Legacy\FP Overplot")

from orbitdata_loading_functions import *
from FP_plot_functions import *
from FP_diagnostics import *
from imagetime_methods import *
from conversion_routines import *
from io_methods import *
from simulation_setup import *
from particle_sim import *

#%%**********************
#FIRST CELL - GET DATA IN
#************************

#choosing comet data to use
inputfilefolder = r"C:\Users\Ollie\Documents\PhD\Comet_Data\Input_files\*pt1.txt"
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
                      horiztag, opts = 'obs')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag,)
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'd10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

#choosing whether to use last used file
last_savename = 'last_image.pickle'
last_savefile = os.path.join(pysav,last_savename)
last_saveexists = os.path.isfile(last_savefile)
getnewfile = True
if last_saveexists == True:
    with open(last_savefile,'rb') as f:
        fitsin = pickle.load(f)
    fitsinfile = os.path.basename(fitsin)
    filebase = fitsinfile[:fitsinfile.find('.')]
    getnewfile = not easygui.ynbox('Use this file: ' + filebase + ' ?')

#choosing fits file to display
if getnewfile == True:
    fitsin = easygui.fileopenbox(default = os.path.join(imagedir,'*'))
    if fitsin == '.': sys.exit("No file selected")

#getting pathnames
fitsinfile = os.path.basename(fitsin)
filebase = fitsinfile[:fitsinfile.find('.')]

#saving name of last used file
with open(last_savefile , 'wb') as f:
    pickle.dump(fitsin, f)

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
    print ("preparing image")
    
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
            if 'A' in obsloc:
                try:
                    w = wcs.WCS(onedimg[0].header, key = 'A')
                except:
                    w = wcs.WCS(onedimg[0].header)
            if 'B' in obsloc:
                w = wcs.WCS(onedimg[0].header)
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
    imgwidth = pixwidth+int(5*border)
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
        
        for xp in range(0,no_points):
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
        
        for xp in range(0,no_points):
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
    
    for div in range(0, (np.size(axisdata[1]))): #RA axis major ticks
        b = d.line([(axisdata[1][div],rdaxis-majt),(axisdata[1][div],rdaxis)],\
        fill = featur_fill)
        tick = str(axisdata[0][div]%360)
        d.text((axisdata[1][div] - len(tick)*5,rdaxis + 10), \
        tick, font=fnt, fill= featur_fill)
        
    for div in range(0, (np.size(axisdata[2]))): #RA axis minor ticks
        b = d.line([(axisdata[2][div],rdaxis-mint),(axisdata[2][div],rdaxis)],\
        fill= featur_fill)
    
    for div in range(0, (np.size(axisdata[4]))): #DEC axis major ticks
        b = d.line([(border+majt,axisdata[4][div]),(border,axisdata[4][div])],\
        fill= featur_fill)
        tick = str(axisdata[3][div])
        d.text((border - 8 - 15*len(tick),axisdata[4][div] - 10 ), \
        tick, font=fnt, fill=featur_fill)
        
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
    
    pix = comimg.load()
    
#%%**************************
#FOURTH CELL - Get Imgtimehdr
#****************************
    
    #some amateur images lack times, need to sometimes use yudish imgtimehdr
    if obsloc == 'Earth':
        
        timemsg = "Choose image time method"
        timechoices = ["Stereo/Soho Filename","Yudish/Earth Filename","User Entry","Yudish Imagetimeheader"]
        reply = easygui.buttonbox(timemsg, choices=timechoices)
        
        if reply == "Stereo/Soho Filename": [ctime,uncertainty_range_exists] = image_time_stereo(filebase)
        if reply == "Yudish/Earth Filename": [ctime,uncertainty_range_exists] = image_time_filename_yuds(filebase)
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
    com_ra_dec = ecliptic_pos2radec(comveceq[comcel - LT_cor,6:9] - obsveceq[comcel,6:9],fixwrapsbool)
    
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
    
    [ltcomcel, vtraj, vtrajcel] = plot_orbit_ecliptic(comobs,comveceq,obsveceq,axisdata,d,comcel,trajfill,ra_img_lower,
    ra_img_higher,border,pixwidth,rafmin,scale,pixheight,decmin,fnt,featur_fill, com_in_image,fixwrapsbool)      
    
    if vtraj is not None:
        plot_orbit_points(d,vtraj,smallfnt,featur_fill,pixheight,pixwidth,border)    
    
    if com_in_image == True:
     
        #setting RGBA colours of lines
        comsunfill = (0,255,0,255)
        manysunfill = (100,255,200,255)
        
        plot_sunearth_vec_ecliptic(d,comveceq,obsveceq,ltcomcel,axisdata,ra_img_higher,ra_img_lower,
                        border,pixwidth,pixheight,rafmin,decmin,scale,comsunfill,featur_fill,fnt,fixwrapsbool)
        
        #plot_many_sunearth_vec(d,comveceq,obsveceq,ctime,ltcomcel,axisdata,ra_img_higher,ra_img_lower,
        #                border,pixwidth,pixheight,rafmin,decmin,scale,manysunfill,featur_fill,fnt,fixwrapsbool)
        
        if obsloc == 'Earth' and uncertainty_range_exists == True:
            plot_compos_unc(d,idls,vtraj,vtrajcel,border,pixwidth,fnt,trajucfill,featur_fill)
    
    else:
        trajucfill = None;comsunfill = None   
    
    #draw the sun if it's in a good place
    sun_pos = ecliptic_pos2radec(-obsveceq[comcel,6:9],fixwrapsbool)
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
    
    with open(picklesavefile, 'wb') as f:
        pickle.dump([comcel, comcel10, ramax, decmax, ramin, decmin, border,
                     pixheight, pixwidth, scale, ctime, dtmin, ra, dec, colr,
                     colb, colg, trajfill, trajucfill, comsunfill, backgr_fill,
                     imgwidth, imgheight, rafmin, rafmax, rapixl, decpixl,
                     com_ra_dec, com_Path, com_in_image, featur_fill, fitscoords], f)
                    
else:
    print ("Loading parameter data")
    
    with open(picklesavefile,'rb') as f:
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
test_mode = True
while test_mode == True:
    
    comimg = Image.open(imgsav)
    d = ImageDraw.Draw(comimg)
    
    simsavefile = os.path.join(pysav, 'simsetupsavs')
    if not os.path.exists(simsavefile ): os.makedirs(simsavefile)
    simsavefile = os.path.join(simsavefile, obsloc)
    if not os.path.exists(simsavefile ): os.makedirs(simsavefile)
    if "Stereo" in obsloc: simsavefile = os.path.join(simsavefile, inst)
    if "Soho" in obsloc: simsavefile = os.path.join(simsavefile, inst)
    if not os.path.exists(simsavefile ): os.makedirs(simsavefile)
    simsavefile  = os.path.join(simsavefile, filebase + '_simsetup.pickle')
    [betau, betal, bno, simtu, simtl, tno, drawopts, sav_bool, test_mode] = simulation_setup(simsavefile)
    
    [ra_m, rafmin, rafmax, fixwrapsbool] = fixwraps(ra, ramax, ramin)    
    
    #finds maximum possible simulateable time and ensure simtu is bounded by this
    simtmax = np.round(ctime.jd - comveceq10[0,0])
    if (simtu > simtmax):
        simtu = simtmax
        print ("Simulation time exceeds orbit data range, setting to max")
    else:
        simtu = simtu
        
    tvals = np.linspace(simtl, simtu, tno)
    
    #create log distributed beta values, accounting that log(0) is impossible
    if (betal == 0):
        bvals = np.logspace(np.log10(0.001), np.log10(betau), num=(bno-1))
        bvals = np.concatenate((np.array([0]), bvals))
        print ("Beta = 0, logarithmic spacing from Beta = 0.001")
    else:
        bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))

    efinp = obsveceq[comcel,6:9]

    if com_in_image == True: image_case = 1
    elif com_in_image == False:
        sync_intersects = np.zeros(tno,dtype = int)
        
        for tval_test in range(0, tno):	
    
                simt10min = int(round(144*tvals[tval_test]))
                pstart = comveceq10[comcel10-simt10min,6:12]
                sim_lo = part_sim(bvals[0],simt10min,30,3,pstart,efinp,dtmin)
                sim_hi = part_sim(bvals[bno-1],simt10min,30,3,pstart,efinp,dtmin)
                ra_dec_lo_test = ecliptic_pos2radec(sim_lo[1][0:3] -
                obsveceq[int(comcel-10*simt10min+sim_lo[0]),6:9],fixwrapsbool)
                ra_dec_hi_test = ecliptic_pos2radec(sim_hi[1][0:3] -
                obsveceq[int(comcel-10*simt10min+sim_hi[0]),6:9],fixwrapsbool)
                sync_box_path = np.array(
                [[ra_dec_lo_test[0],ra_dec_lo_test[1]],
                [ra_dec_hi_test[0],ra_dec_hi_test[1]]])
                sync_Path = mplPath.Path(sync_box_path)
                sync_intersects[tval_test] = com_Path.intersects_path(sync_Path)

        if (np.sum(sync_intersects) > 0): image_case = 2
        elif (np.sum(sync_intersects) == 0): image_case = 3        
    
    #prep for simulation loop
    simres = np.zeros((tno,bno,18),dtype = float)
    bmax = np.zeros((tno),dtype = int)
    bmin = np.zeros((tno),dtype = int)
    tmax = np.zeros((bno),dtype = int)
    tmin = np.zeros((bno),dtype = int)
    tidx = 0

#%%*********************************
#EIGHT CELL - Simulating Dust Motion
#***********************************
    if image_case == 1:
        
        while (tidx < tno):
            bidx = 0
            point_in_image = 1
            simt10min = int(round(144*tvals[tidx]))
            pstart = comveceq10[comcel10-simt10min,6:12]
            tswap = 2454112.6 - comveceq10[comcel10-simt10min,0]
            while (bidx < bno and point_in_image == 1):
                sim = part_sim_lorentz(bvals[bidx],simt10min,50,1,pstart,efinp,dtmin,tswap)
                simres[tidx,bidx,0] = float(simt10min)/144
                simres[tidx,bidx,1] = bvals[bidx]
                simres[tidx,bidx,2] = sim[0] #length of simulation in minutes
                simres[tidx,bidx,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
                simres[tidx,bidx,4:10] = sim[1] #finishing pos/vel
                simres[tidx,bidx,10:12] = ecliptic_pos2radec(sim[1][0:3] - 
                obsveceq[int(simres[tidx,bidx,3]),6:9],fixwrapsbool)
                point_in_image = int(com_Path.contains_point(
                (simres[tidx,bidx,10],simres[tidx,bidx,11])))
                simres[tidx,bidx,12] = ra2xpix(simres[tidx,bidx,10],border,pixwidth,rafmin,scale)
                simres[tidx,bidx,13] = dec2ypix(simres[tidx,bidx,11],border,pixheight,decmin,scale)                 
                simres[tidx,bidx,14] = point_in_image
                simres[tidx,bidx,17] = getdustphase(sim[1][0:3],obsveceq[int(simres[tidx,bidx,3]),6:9])
                bidx += 1
            try:simres[tidx,bidx - 1 + point_in_image,15] = 0
            except: pass
            try:bmax[tidx] = simres[tidx,:,14].nonzero()[0][-1]
            except: bmax[tidx] = 0            
            tidx += 1
            print (float(tidx)*100/tno)
            
        simres[np.where(simres[:,:,14]==0)] = 0
        bidx_list = np.arange(bno)
        deletions = 0
        for bidx in range(0,bno):
            [minblock,maxblock] = find_largest_nonzero_block(simres[:,bidx,14])
            if minblock == None:
                bidx_list = np.delete(bidx_list,bidx - deletions)
                deletions += 1
            else:
                tmin[bidx] = minblock
                tmax[bidx] = maxblock        
                
        tidx_list = np.unique(np.nonzero(simres[:,:,14])[0])
        
        tidx_tt = np.zeros_like(tvals,dtype = int)        
        bidx_tt = np.zeros_like(bvals,dtype = int) 
        tidx_tt[tidx_list] = 1
        bidx_tt[bidx_list] = 1
        
    if image_case == 2:
        
        tidx_list = np.arange(np.where(sync_intersects == 1)[0][0],
                              np.where(sync_intersects == 1)[0][-1]+1)
        syn_exit_tt = np.zeros((2,2),dtype = int)
        syn_exit_tt[0,1] = 1
        bprev = bno-1
        
        for tidx_list_val in range(0,np.size(tidx_list)):
            bidx = bno-1
            tidx = tidx_list[tidx_list_val]
            syn_exited_image = 0
            prev_stat = 0
            simt10min = int(round(144*tvals[tidx]))  
            pstart = comveceq10[comcel10-simt10min,6:12]
            tswap = 2454112.6 - comveceq10[comcel10-simt10min,0]
            while (bidx >= 0 and (syn_exited_image == 0 or bidx >= bprev)):
                sim = part_sim_lorentz(bvals[bidx],simt10min,50,3,pstart,efinp,dtmin,tswap)
                simres[tidx,bidx,0] = float(simt10min)/144
                simres[tidx,bidx,1] = bvals[bidx]
                simres[tidx,bidx,2] = sim[0] #length of simulation in minutes
                simres[tidx,bidx,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
                simres[tidx,bidx,4:10] = sim[1] #finishing pos/vel
                simres[tidx,bidx,10:12] = ecliptic_pos2radec(sim[1][0:3] - 
                obsveceq[int(simres[tidx,bidx,3]),6:9],fixwrapsbool)
                point_in_image = int(com_Path.contains_point
                ((simres[tidx,bidx,10],simres[tidx,bidx,11])))
                simres[tidx,bidx,12] = ra2xpix(simres[tidx,bidx,10],border,pixwidth,rafmin,scale)
                simres[tidx,bidx,13] = dec2ypix(simres[tidx,bidx,11],border,pixheight,decmin,scale)                 
                simres[tidx,bidx,14] = point_in_image
                syn_exited_image = syn_exit_tt[point_in_image,prev_stat]
                prev_stat = point_in_image
                simres[tidx,bidx,17] = getdustphase(sim[1][0:3],obsveceq[int(simres[tidx,bidx,3]),6:9])
                bidx -= 1
            print (float(tidx)*100/tno)
            bprev = bidx
            try:simres[tidx,bidx+1-point_in_image,15] = 0
            except:pass
            try:
                bmin[tidx] = simres[tidx,:,14].nonzero()[0][0]
                bmax[tidx] = simres[tidx,:,14].nonzero()[0][-1]
            except: sync_intersects[tidx] = 0
        
        simres[np.where(simres[:,:,14]==0)] = 0
        bidx_list = np.arange(bno)
        deletions = 0
        for bidx in range(0,bno):
            [minblock,maxblock] = find_largest_nonzero_block(simres[:,bidx,14])
            if minblock == None:
                bidx_list = np.delete(bidx_list,bidx - deletions)
                deletions += 1
            else:
                tmin[bidx] = minblock
                tmax[bidx] = maxblock

        tidx_tt = np.zeros_like(tvals,dtype = int)        
        bidx_tt = np.zeros_like(bvals,dtype = int) 
        tidx_tt[tidx_list] = 1
        bidx_tt[bidx_list] = 1         
                        
    if simres[:,:,12].max() == 0:
        sys.exit("Image " + fitsinfile + " was no good.")
        image_case = 3
            
    if image_case == 3:
        sys.exit("Image " + fitsinfile + " was no good.")
        
    non_zeros_0 = simres[:,:,14].nonzero()[0]
    non_zeros_1 = simres[:,:,14].nonzero()[1]
    no_points = non_zeros_0.size
    nu_radecs = np.zeros((no_points,2))
    
    for x in range(0,no_points):
        nu_radecs[x,0] = simres[non_zeros_0[x],non_zeros_1[x],10]
        nu_radecs[x,1] = simres[non_zeros_0[x],non_zeros_1[x],11]
        
    if np.mean(ra_m-ra) == 360:
        nu_radecs[:,0] = nu_radecs[:,0] - 360
   
    try:
        pixlocs = w.wcs_world2pix(nu_radecs,0)
    except:
        onedimg = fits.open(fitscoords)
        w = wcs.WCS(onedimg[0].header)
        pixlocs = w.wcs_world2pix(nu_radecs,0)
    
    for x in range(0,no_points):
        simres[non_zeros_0[x],non_zeros_1[x],15] = pixlocs[x][0]
        simres[non_zeros_0[x],non_zeros_1[x],16] = pixlocs[x][1]
        
    if (sav_bool == True):
        simressavefile = os.path.join(pysav, 'simres')
        if not os.path.exists(simressavefile): os.makedirs(simressavefile)
        simressavefile = os.path.join(simressavefile, obsloc)
        if not os.path.exists(simressavefile): os.makedirs(simressavefile)
        if "Stereo_A" in obsloc: simsavbase = filebase[:filebase.find('A')+1]
        elif "Stereo_B" in obsloc: simsavbase = filebase[:filebase.find('B')+1]
        elif "Soho" in obsloc: simsavbase = filebase[:filebase.find('Clear')+5]
        else: simsavbase = filebase
        simressavefile = os.path.join(simressavefile, simsavbase + '_' + str(betal) + '_'
                         + str(betau) + '_' + str(bno)+ '_' + str(simtl) + '_'
                         + str(simtu) + '_' + str(tno))
        simressavefile = simressavefile.replace('.','\'')
        np.save(simressavefile, simres)
        with open(simressavefile + '_parameters.pickle' , 'wb') as f:
            pickle.dump([tmax, bmax, tvals, bvals], f)
    
#%%***************************
#NINTH CELL - Plot dust motion
#*****************************

    dynfill = (255,0,0,255) #syndynes
    chrfill = (255,192,0,255) #synchrones
    drfill = (255,0,255,255) #data points and data
    
    fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
    fnt = ImageFont.truetype(fontloc, 20)  
    bt_anno_idx = 0
    
    if "Wide" in drawopts:
        tspacing = int(tno/5); bspacing = int(bno/5)
    else:
        tspacing = 1; bspacing = 1
    if "Syndynes" in drawopts: draw_syndynes(dynfill,d,simres,bno,rapixl,decpixl,tmin,tmax,bidx_list,bspacing,bvals)
    if "Synchrones" in drawopts: draw_synchrones(chrfill,d,simres,tno,rapixl,decpixl,bmin,bmax,tidx_list,tspacing,tvals)
    if "Data Points" in drawopts: draw_datap(drfill,d,simres)
    if "Dust Phase Angles" in drawopts: [smax, smin, grad, colormap] = draw_phase_points(d,simres)
    if "Data Region Enclosed" in drawopts: draw_data_reg(drfill,d,simres,bmax,bmin,bidx_list,tmax,tmin,tidx_list,border,pixwidth)

    if not "Dust Phase Angles" in drawopts:
        bt_anno_idx = annotate_plotting(d,drawopts,border,pixwidth,fnt,featur_fill,dynfill,chrfill,drfill)        
    else:
        bt_anno_idx = annotate_dustphase(d,border,pixwidth,featur_fill,fnt,smax,smin,grad,colormap)
        
    write_bt_ranges(d,border, pixwidth, fnt, featur_fill,
                    betau, betal, simtu, simtl, bt_anno_idx)

    if (drawopts != "No Image"):
        cimgdir = os.path.join(imagedir, 'FPplots')
        if not os.path.exists(cimgdir): os.makedirs(cimgdir)
        cimgsav = os.path.join(cimgdir, 'FP_lorentz_' + filebase + '.png')    
        comimg.save(cimgsav,'png')
        webbrowser.open(cimgsav)