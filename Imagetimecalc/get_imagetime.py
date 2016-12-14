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

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from orbitdata_loading_functions import orb_vector, orb_obs
from FP_plot_functions import ra2xpix, dec2ypix, setaxisup, plotpixel
from conversion_routines import pos2radec, fixwraps, find_largest_nonzero_block 
from io_methods import correct_for_imagetype

#%%**********************
#FIRST CELL - GET DATA IN
#************************

if False:
    #choosing comet data to use
    inputfilefolder = "C:\PhD\Comet_data\Input_files\*pt1.txt"
    inputfile = easygui.fileopenbox(default = inputfilefolder)
    
    #reading main comet parameters
    with open(inputfile, "r") as c:
        cdata = c.readlines()
        comname = cdata[30][12:]
        comdenom = cdata[31][13:-2]
        orbitdir = cdata[25][23:-2]
        timedir = cdata[28][26:-2]
        pysav = cdata[27][24:-2]
        horiztag = cdata[40][10:]
    
    obsloc = 'Earth'    
    #import the orbit data
    obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                          horiztag, opts = 'obs,eq')
    comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                          horiztag, opts = 'eq')
    comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                            horiztag, opts = 'eq,d10')
    comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

#choosing fits file to display and getting pathnames
fitsin = easygui.fileopenbox(default = os.path.join(timedir,'*'))
if fitsin == '.': sys.exit("No file selected")
fitsinfile = os.path.basename(fitsin)
filebase = fitsinfile[:string.find(fitsinfile,'.')]

#ensures image inputted correctly depending on size of data cube
[colr, colg, colb, fitscoords] = correct_for_imagetype(timedir, fitsin, fitsinfile)

        
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
radecs = w.wcs_pix2world(coords, 0)
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
comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,(0,0,0,255))
d = ImageDraw.Draw(comimg)
for x in xrange(0,np.shape(colr)[0]-1):
    for y in xrange (0,np.shape(colr)[1]-1):
        plotpixel(d,x,y,ra_m,dec,border,pixwidth,pixheight,decmin,rafmin,
                  scale,colr[x,y],colg[x,y],colb[x,y])
            
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

comlocra = ra2xpix(comobs[:,5], border, pixwidth, rafmin, scale)
comlocdec = dec2ypix(comobs[:,6], border, pixheight, decmin, scale)

decgoodlocs = np.where((comlocdec < (imgheight-border*1.5)) & (comlocdec > border*1.5))[0]
ragoodlocs = np.where((comlocra < (imgwidth-border*2.5)) & (comlocra > border*1.5))[0]
goodlocs = np.intersect1d(decgoodlocs,ragoodlocs)

for o in xrange(0,np.size(goodlocs)-1):
    ocel = goodlocs[o]
    d.line( [ (comlocra[ocel] ,comlocdec[ocel]), 
    (comlocra[ocel+1] ,comlocdec[ocel+1]) ],
    fill = (255,0,0,255))
    
imgsave = os.path.join(os.path.join(timedir,'images'),
                       filebase.split('_')[-1] + '_astrometrytest.png')
comimg.save(imgsave,'png')

print "Done!"