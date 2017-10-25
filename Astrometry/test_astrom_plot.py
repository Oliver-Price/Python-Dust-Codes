import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont
import datetime
import astropy.time
import easygui

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from io_methods import get_obs_loc, get_stereo_instrument, get_soho_instrument, get_hih_low
from orbitdata_loading_functions import orb_vector, orb_obs
from conversion_routines import fixwraps, pos2radec
from FP_plot_functions import setaxisup, plotpixel2, ra2xpix, dec2ypix
from imagetime_methods import image_time_stereo

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

#rdcsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\rdc.npy'

outdir = os.path.join(imagedir,'aligned_pngs')
if not os.path.exists(outdir): os.makedirs(outdir)

dir_list = sorted(os.listdir(imagedir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)
#rc = np.load(rdcsav)

if "C2006P1" in imagedir:
    if obsloc == "Soho": 
        trafmin = 285.493251302; trafmax = 306.277182607
        tdecmin = -29.4102183696; tdecmax = -12.1504198155
    if obsloc == "Stereo_A":
        if inst == 'HI-1':
            trafmax = 332.57662347125824; trafmin = 290.18476496395908
            tdecmax = 2.5007098119796005; tdecmin = -28.041542198531385
        if inst == "HI-2":
            trafmax = 404.97809546316046; trafmin = 289.36252667052929
            tdecmax = 58.864096027413353; tdecmin = -37.364638404644303
    if obsloc == "Stereo_B":
        if inst == "HI-1":
            trafmax = 404.97809546316046; trafmin = 282.46071535529717
            tdecmax = 58.864096027413353; tdecmin = -37.364638404644303
        if inst == "HI-2":
            trafmax = 404.97809546316046; trafmin = 224.39349694174101
            tdecmax = 69.369679615559832; tdecmin = -37.364638404644303
    
elif "C2011L4" in imagedir:
    trafmin = 180.35639546709564; trafmax = 235.62536933813792
    tdecmin = -38.206859795611066; tdecmax = 13.843362244035761
    
elif "C2002V1" in imagedir:
    trafmin = 318.99291080194598; trafmax = 344.28481924562749
    tdecmin = -22.783497571686073; tdecmax = -0.99477817124254442
    
elif "C2011W3" in imagedir:
    trafmin = 252.62004219486914; trafmax = 273.42781500602257
    tdecmin = -31.034491925907702; tdecmax = -14.431726967874498

for fits_no in range(0,fits_total):

    print (fits_no)
    image_basename = fits_list[fits_no]
    fitsin = os.path.join(imagedir, fits_list[fits_no])
    imgsav = os.path.join(outdir, image_basename.split('.')[0] + '.png')
    if not os.path.exists(imgsav):
        
        onedimg = fits.open(fitsin)
        try:
            w = wcs.WCS(onedimg[0].header, key = 'A', naxis=2)
        except:
            w = wcs.WCS(onedimg[0].header, naxis=2)
        colours = np.nan_to_num(onedimg[0].data)
        
        #make a 2xN array of all pixel locations
        ya = colours.shape[0]
        xa = colours.shape[1]
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
        
        [ra_m, rafmin, rafmax, bool_val] = fixwraps(ra, ramax, ramin)
        
        backgr_fill = (0,0,0,255)    
        featur_fill = (255,255,255,255)
        
        pixheight = 800
        pixwidth = int(pixheight*(trafmax - trafmin)/(tdecmax - tdecmin))
        border = 100
        scale = pixheight/(tdecmax - tdecmin)
        imgwidth = pixwidth+int(4*border)
        imgheight = pixheight+int(3*border)
        comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,backgr_fill)
        d = ImageDraw.Draw(comimg)
        
    #    for rdc in xrange(0, fits_total):
    #        r1 = ra2xpix(rc[rdc,0], border, pixwidth, trafmin, scale)
    #        r2 = ra2xpix(rc[rdc,1], border, pixwidth, trafmin, scale)
    #        r3 = ra2xpix(rc[rdc,2], border, pixwidth, trafmin, scale)
    #        r4 = ra2xpix(rc[rdc,3], border, pixwidth, trafmin, scale)
    #        d1 = dec2ypix(rc[rdc,4], border, pixheight, tdecmin, scale)
    #        d2 = dec2ypix(rc[rdc,5], border, pixheight, tdecmin, scale)
    #        d3 = dec2ypix(rc[rdc,6], border, pixheight, tdecmin, scale)
    #        d4 = dec2ypix(rc[rdc,7], border, pixheight, tdecmin, scale)
    #        d.point([(r1,d1),(r2,d2),(r3,d3),(r4,d4)],(255,0,0,255))
          
    
        [hih,low] = get_hih_low(comdenom,obsloc,inst)
        
        grad = (255 / (hih - low)    )  
        
        colcr = (np.clip(np.round((np.clip(colours,1e-20,999999999)- low)*grad),0,255)).astype(int)
        colcg = (np.clip(np.round((np.clip(colours,1e-20,999999999)- low)*grad),0,255)).astype(int)
        colcb = (np.clip(np.round((np.clip(colours,1e-20,999999999)- low)*grad),0,255)).astype(int)
        
        imagemask = np.ones_like(colours)[:-1,:-1]
        
        non_zeros_0 = np.where(imagemask == 1)[0]
        non_zeros_1 = np.where(imagemask == 1)[1]
        no_points = non_zeros_0.size 
        
        rp = ra2xpix(ra_m, border, pixwidth, trafmin, scale)
        dp = dec2ypix(dec, border, pixheight, tdecmin, scale)
        
        for xp in range(0,no_points):
            x = non_zeros_0[xp]
            y = non_zeros_1[xp]
            plotpixel2(d,x,y,rp,dp,colcr[x,y],colcg[x,y],colcb[x,y])
        
        if 'S' in obsloc:
            [ctime,uncertainty_range_exists] = image_time_stereo(image_basename)
        
        comcel = np.where(abs(obsveceq[:,0] - ctime.jd)==abs(obsveceq[:,0] - ctime.jd).min())[0][0]
    
        #find rough cell range of traj from observer data
        obsraloc = np.where((comobs[:,5] > trafmin) \
        & (comobs[:,5] < trafmax))[0]
        obsdecloc = np.where((comobs[:,6] > tdecmin) \
        & (comobs[:,6] < tdecmax))[0]
        trajrough = np.intersect1d(obsraloc,obsdecloc)
        
        if trajrough.size != 0:    
        
            #use this to calculate ra and dec of comet for a purposely oversized range
            vno = 0; vext = int(np.size(trajrough))
            vtraj = np.empty((np.size(trajrough)+2*vext-1,11),dtype = float)
            tcellmax = min(trajrough[-1] + vext, np.shape(comveceq)[0])
            for tcell in range(trajrough[0] - vext,tcellmax):
                vtemp = comveceq[tcell,6:9] - obsveceq[comcel,6:9]    
                ptemp = pos2radec(vtemp,bool_val)
                vtraj[vno,0] = ptemp[0]
                vtraj[vno,1] = ptemp[1]
                vtraj[vno,4] = tcell
                vtraj[vno,5] = tcell + round(np.linalg.norm(vtemp)*8.316746397269274)
                vtraj[vno,6:11] = comveceq[tcell,1:6]
                vno +=1
            
        #use these ra and dec values to slim data down to within image borders
        trajrange = np.intersect1d(
                    np.intersect1d( np.where(vtraj[:,0] < trafmax)[0],
                                    np.where(vtraj[:,0] > trafmin)[0]),
                    np.intersect1d( np.where(vtraj[:,1] < tdecmax)[0],
                                    np.where(vtraj[:,1] > tdecmin)[0]))                  
        vtraj = vtraj[trajrange[0]:trajrange[-1],:]
        
        #convert to ra and dec, and plot
        vtraj[:,2] = ra2xpix(vtraj[:,0],border,pixwidth,trafmin,scale)
        vtraj[:,3] = dec2ypix(vtraj[:,1],border,pixheight,tdecmin,scale)
        for ta in range(0, (np.shape(vtraj)[0]-1)):
            d.line([(vtraj[ta,2],vtraj[ta,3]),(vtraj[ta+1,2],vtraj[ta+1,3])],\
            fill = (0,255,255,255))
        
        #draws a border       
        a = d.polygon([(border,border),(border*2+pixwidth,border), \
        (border*2+pixwidth,border*2+pixheight),(border,border*2+pixheight)], \
        outline = featur_fill)
        
        #most of the dirty stuff is bunged into this function
        axisdata = setaxisup(trafmax,trafmin,tdecmax,tdecmin,border,pixheight,pixwidth,scale)
    
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
        
        plttitle = ('Soho Lasco C3 Clear at: ' + ctime.isot[0:16].replace('T',' at '))
        d.text((1.5*border + pixwidth*0.5 - len(plttitle)*5 - 60,.35*border), \
        plttitle, font=largefnt, fill= featur_fill)
        
        comimg.save(imgsav,'png')
