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
import astropy
import datetime
import matplotlib.path as mplPath
from PIL import Image, ImageDraw, ImageFont

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

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
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)
venobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag, venusmode = True)
mercobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag, mercurymode = True)

dir_list = sorted(os.listdir(imagedir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

if not os.path.exists(os.path.join(imagedir,'imgtime_test_venmerc')): os.makedirs(os.path.join(imagedir,'imgtime_test_venmerc'))
        
for fidx in range (100,101):#fits_total):
    
    fitsinfile = fits_list[fidx]
    filebase = fits_list[fidx].split('.')[0]
    fitsin = os.path.join(imagedir,fitsinfile)
    print (filebase)
        
    imgsave = os.path.join(os.path.join(imagedir,'imgtime_test_venmerc'),
                           filebase + '_astrometrytest.png')
    
    if not os.path.exists(imgsave):
           
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
                try:
                    w = wcs.WCS(onedimg[0].header, key = 'A')
                except:
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
        
        [hih,low] = get_hih_low(comdenom,obsloc,inst)
        
    #    hih = 325000
    #    low = 40000
        
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
        
        #comet
        comlocra = ra2xpix(comobs[:,5], border, pixwidth, rafmin, scale)
        comlocdec = dec2ypix(comobs[:,6], border, pixheight, decmin, scale)
        
        decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
        ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
        goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)
        
        for o in range(0,np.size(goodlocs)-1):
            ocel = goodlocs[o]
            d.line( [ (comlocra[ocel] ,comlocdec[ocel]), 
            (comlocra[ocel+1] ,comlocdec[ocel+1]) ],
            fill = (255,0,0,255))
        
        #venus
        comlocra = ra2xpix(venobs[:,5], border, pixwidth, rafmin, scale)
        comlocdec = dec2ypix(venobs[:,6], border, pixheight, decmin, scale)
        
        decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
        ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
        goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)
        
        for o in range(0,np.size(goodlocs)-1):
            ocel = goodlocs[o]
            d.line( [ (comlocra[ocel] ,comlocdec[ocel]), 
            (comlocra[ocel+1] ,comlocdec[ocel+1]) ],
            fill = (0,255,0,255))
        
        #merc
        comlocra = ra2xpix(mercobs[:,5], border, pixwidth, rafmin, scale)
        comlocdec = dec2ypix(mercobs[:,6], border, pixheight, decmin, scale)
        
        decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
        ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
        goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)
        
        for o in range(0,np.size(goodlocs)-1):
            ocel = goodlocs[o]
            d.line( [ (comlocra[ocel] ,comlocdec[ocel]), 
            (comlocra[ocel+1] ,comlocdec[ocel+1]) ],
            fill = (0,0,255,255))
        
        #imagetime  
        csec = int(filebase[13:15])
        cmin = int(filebase[11:13])
        chour = int(filebase[9:11])
        cday = int(filebase[6:8])
        cmonth = int(filebase[4:6])
        cyear = int(filebase[0:4])
        ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                    chour , cmin, csec))
        
        orbit_cells = np.where(venobs[:,1] == cmonth)[0]
        orbit_cells = np.intersect1d(np.where(venobs[:,2] == cday)[0],orbit_cells)
        orbit_cells = np.intersect1d(np.where(venobs[:,3] == chour)[0],orbit_cells)
        comcell = np.intersect1d(np.where(venobs[:,4] == cmin)[0],orbit_cells)[0]
        
        xsiz = 5
        com_ra_loc = ra2xpix(comobs[comcell,5], border, pixwidth, rafmin, scale)
        com_dec_loc = dec2ypix(comobs[comcell,6], border, pixheight, decmin, scale)
        d.line( [ ( com_ra_loc - xsiz , com_dec_loc - xsiz ) ,
                  ( com_ra_loc + xsiz , com_dec_loc + xsiz ) ] ,
                  fill = (255,0,255,255) )  
        d.line( [ ( com_ra_loc - xsiz , com_dec_loc + xsiz ) ,
                  ( com_ra_loc + xsiz , com_dec_loc - xsiz ) ] ,
                  fill = (255,0,255,255) )
    
        com_ra_loc = ra2xpix(venobs[comcell,5], border, pixwidth, rafmin, scale)
        com_dec_loc = dec2ypix(venobs[comcell,6], border, pixheight, decmin, scale)
        d.line( [ ( com_ra_loc - xsiz , com_dec_loc - xsiz ) ,
                  ( com_ra_loc + xsiz , com_dec_loc + xsiz ) ] ,
                  fill = (0,255,255,255) )  
        d.line( [ ( com_ra_loc - xsiz , com_dec_loc + xsiz ) ,
                  ( com_ra_loc + xsiz , com_dec_loc - xsiz ) ] ,
                  fill = (0,255,255,255) )
        
        com_ra_loc = ra2xpix(mercobs[comcell,5], border, pixwidth, rafmin, scale)
        com_dec_loc = dec2ypix(mercobs[comcell,6], border, pixheight, decmin, scale)
        d.line( [ ( com_ra_loc - xsiz , com_dec_loc - xsiz ) ,
                  ( com_ra_loc + xsiz , com_dec_loc + xsiz ) ] ,
                  fill = (255,255,0,255) )  
        d.line( [ ( com_ra_loc - xsiz , com_dec_loc + xsiz ) ,
                  ( com_ra_loc + xsiz , com_dec_loc - xsiz ) ] ,
                  fill = (255,255,0,255) )
        
        comimg.save(imgsave,'png')
