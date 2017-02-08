import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont
import datetime
import astropy.time
import time

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from conversion_routines import fixwraps
from FP_plot_functions import setaxisup, plotpixel2, ra2xpix, dec2ypix

fitsdir = r'C:\PhD\Comet_data\Comet_NEAT_C2002V1\Gallery\Soho\C3_Clear'
outdir = r'C:\PhD\Comet_data\Comet_NEAT_C2002V1\Gallery\Soho\C3_Clear\aligned_pngs'
rdcsav = r'C:\PhD\Comet_data\Comet_NEAT_C2002V1\Gallery\Soho\rdc.npy'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".fits" in s]
fits_total = len(fits_list)
rc = np.load(rdcsav)

for fits_no in xrange(0,fits_total):

    print fits_no
    image_basename = fits_list[fits_no]
    fitsin = os.path.join(fitsdir, fits_list[fits_no])
    imgsav = os.path.join(outdir, image_basename.split('.')[0] + '.png')
    if not os.path.exists(imgsav):
        onedimg = fits.open(fitsin)
        w = wcs.WCS(onedimg[0].header, naxis=2)
        colours = onedimg[0].data
        
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
        
        if "C2006P1" in fitsdir:
            trafmin = 285.493251302; trafmax = 306.277182607
            tdecmin = -29.4102183696; tdecmax = -12.1504198155
            
        if "C2011L4" in fitsdir:
            trafmin = 180.35639546709564; trafmax = 235.62536933813792
            tdecmin = -38.206859795611066; tdecmax = 13.843362244035761
            
        if "C2002V1" in fitsdir:
            trafmin = 318.99291080194598; trafmax = 344.28481924562749
            tdecmin = -22.783497571686073; tdecmax = -0.99477817124254442 
        
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
          
    
        low = 1e-13; hih = 1e-11
        
        grad = (255 / (np.log10(hih) - np.log10(low)))      
        
        colcr = (np.clip(np.round((np.log10(np.clip(colours,1e-20,999999999))- np.log10(low))*grad),0,255)).astype(int)
        colcg = (np.clip(np.round((np.log10(np.clip(colours,1e-20,999999999))- np.log10(low))*grad),0,255)).astype(int)
        colcb = (np.clip(np.round((np.log10(np.clip(colours,1e-20,999999999))- np.log10(low))*grad),0,255)).astype(int)
        
        imagemask = np.ones_like(colours)[:-1,:-1]
        
        non_zeros_0 = np.where(imagemask == 1)[0]
        non_zeros_1 = np.where(imagemask == 1)[1]
        no_points = non_zeros_0.size 
        
        rp = ra2xpix(ra_m, border, pixwidth, trafmin, scale)
        dp = dec2ypix(dec, border, pixheight, tdecmin, scale)
        
        for xp in xrange(0,no_points):
            x = non_zeros_0[xp]
            y = non_zeros_1[xp]
            plotpixel2(d,x,y,rp,dp,colcr[x,y],colcg[x,y],colcb[x,y])
                    
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
            
        csec = int(image_basename[13:15])
        cmin = int(image_basename[11:13])
        chour = int(image_basename[9:11])
        cday = int(image_basename[6:8])
        cmonth = int(image_basename[4:6])
        cyear = int(image_basename[0:4])
        ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday, chour , cmin, csec))
        
        #axis labels
        d.text((1.5*border + pixwidth*0.5 - 145,pixheight + 2.5*border), \
        "Right Ascension (Degrees)", font=fnt, fill= featur_fill)
        d.text((0.25*border - 10,0.75*border - 20), \
        "Declination (Degrees)", font=fnt, fill= featur_fill)
        
        plttitle = ('Soho Lasco C3 Clear at: ' + ctime.isot[0:16].replace('T',' at '))
        d.text((1.5*border + pixwidth*0.5 - len(plttitle)*5 - 60,.35*border), \
        plttitle, font=largefnt, fill= featur_fill)
        
        comimg.save(imgsav,'png')
