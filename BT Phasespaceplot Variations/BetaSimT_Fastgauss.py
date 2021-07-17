#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************

import easygui
import os
import sys
import pickle
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont
import webbrowser
import datetime
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")  

from BT_io_functions import fixed_image_times,dategetter
from BT_plot_functions import plotpixel, logplotpixel, bt_setaxisup
from io_methods import get_obs_loc, get_stereo_instrument, get_soho_instrument, get_hih_low
from conversion_routines import fixwraps

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
    
#choosing fits file to display and getting pathnames
pngdir = os.path.join(imagedir, 'cometplots')
pngin = easygui.fileopenbox(default = os.path.join(pngdir,'*'))
infile = os.path.basename(pngin)
filebase = infile[:infile.find('.')]
if 'inverted' in infile:
    filebase = filebase.split("_inverted")[0]

#parameter savefile location
picklesavefile = os.path.join(pysav, 'imgsavs')
picklesavefile = os.path.join(picklesavefile, obsloc)
if "S" in obsloc: picklesavefile = os.path.join(picklesavefile, inst)
picklesavefile = os.path.join(picklesavefile, filebase + '_plot_param.pickle')

greyscale = True

#check vital information exists
if not os.path.exists(picklesavefile):
    sys.exit("Image not calculated")

#import important information
with open(picklesavefile,'rb') as f:
        dparameters = pickle.load(f)
        comcel = dparameters[0]
        comcel10 = dparameters[1]
        ramax = dparameters[2]
        decmax = dparameters[3]
        ramin = dparameters[4]
        decmin = dparameters[5]        
        ctime = dparameters[10]
        dtmin = dparameters[11]
        ra = dparameters[12]
        dec = dparameters[13]
        colr = dparameters[14]
        colb = dparameters[15]
        colg = dparameters[16]
        rapixl = dparameters[25]
        decpixl = dparameters[26]
        com_ra_dec = dparameters[27]
        fitscoords = dparameters[31]
        
#find and import simulation results
simresdir = os.path.join(pysav, 'simres')
if "Stereo" in obsloc: simsavbase = filebase[:filebase.find('A')+1]
elif "Soho" in obsloc:
    if inst == 'C3_Clear':
        simsavbase = filebase[:filebase.find('Clear')+5]
    if inst == 'C3_Blue':
        simsavbase = filebase[:filebase.find('Bl')+2]   
else: simsavbase = filebase
simresdir = os.path.join(simresdir, obsloc)
simin = easygui.fileopenbox(default = os.path.join(simresdir, simsavbase + '*'))
simres = np.load(simin)

with open(simin[:-4] + '_parameters.pickle','rb') as f:
    sparameters = pickle.load(f)
    tmax = sparameters[0]
    bmax = sparameters[1]
    tvals = sparameters[2]
    bvals = sparameters[3]

#%%*********************************************************
#SECOND CELL - IDENTIFYING PIXEL VALUES TO USE FOR BETA/SIMT
#***********************************************************

onedimg = fits.open(fitscoords)
if 'Earth' in obsloc:
    w = wcs.WCS(onedimg[0].header)
elif 'Stereo' in obsloc:
    if comdenom == 'c2011l4':
        w = wcs.WCS(onedimg[0].header)
    elif comdenom == 'c2006p1':
        w = wcs.WCS(onedimg[0].header, key = 'A')
elif 'Soho' in obsloc:
    w = wcs.WCS(onedimg[0].header)

#putting ra into sensible values
[ra_m, rafmin, rafmax, fixwrapsbool] = fixwraps(ra, ramax, ramin)
    
#imagemask for HI-2
if "Stereo" in obsloc:
    if '2' in inst:
        maskedimg = fits.open('C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\hi2A_mask.fts')
        imagemask  = maskedimg[0].data[::2,::2]
        
    elif '1' in inst:
        direndindx = os.path.dirname(fitscoords).find("HI-1") + 4
        fits_dir = os.path.dirname(fitscoords)[:direndindx]
        fits_list = os.listdir(fits_dir)
        fitstr = [s for s in fits_list if filebase[:21] in s][0]
        orig_fits_file = os.path.join(fits_dir,fitstr)
        
        hdulist2 = fits.open(orig_fits_file)
        colours = (hdulist2[0].data)
        
        imagemask = np.ones_like(colr)
        imagemask[np.where(colours > 1.5e6)] = 0
        
elif "Soho" in obsloc:
        if "Blue" in inst:
            direndindx = os.path.dirname(fitscoords).find("Blue") + 4
            fits_dir = os.path.dirname(fitscoords)[:direndindx]
            fits_list = os.listdir(fits_dir)
            fitstr = [s for s in fits_list if filebase[:18] in s][0]
            orig_fits_file = os.path.join(fits_dir,fitstr)
        elif "Clear" in inst:
            direndindx = os.path.dirname(fitscoords).find("Clear") + 5
            fits_dir = os.path.dirname(fitscoords)[:direndindx]
            fits_list = os.listdir(fits_dir)
            fitstr = [s for s in fits_list if filebase[:21] in s][0]
            orig_fits_file = os.path.join(fits_dir,fitstr)
        
        hdulist2 = fits.open(orig_fits_file)
        colours = (hdulist2[0].data)
        imagemask = np.ones_like(colr)
        #imagemask[np.where(colours > 6000)] = 0 #1300 #6000
        
else: imagemask = np.ones_like(colr)

#check if this has already been done
colmapsav = simin[:-4] + '_srcolors.npy'
locmapsav = simin[:-4] + '_pixelmapping.txt'

tno = np.size(tvals); bno = np.size(bvals)
srcolors = np.zeros((tno,bno,6),dtype=float)
srcolors_old = np.zeros((tno,bno,4),dtype=float)

simres_floor = np.floor(simres[:,:,14:17]).astype(int)
simres_rounded = np.clip(np.round(simres[:,:,15:17]).astype(int), 0, 5000)

non_zeros_0 = simres_floor[:,:,0].nonzero()[0]
non_zeros_1 = simres_floor[:,:,0].nonzero()[1]
no_points = non_zeros_0.size

for x in range(0,no_points):
    tidx = non_zeros_0[x]
    bidx = non_zeros_1[x]
    dval = simres[tidx,bidx,15]
    rval = simres[tidx,bidx,16]
    drou = simres_rounded[tidx,bidx,0]
    rrou = simres_rounded[tidx,bidx,1]
    dflo = simres_floor[tidx,bidx,1]
    rflo = simres_floor[tidx,bidx,2]
    srcolors[tidx,bidx,0] = imagemask[rrou,drou]
    srcolors[tidx,bidx,1] = (colr[rflo,dflo]*(rflo - rval + 1)*(dflo - dval + 1) + 
                             colr[rflo+1,dflo]*(rval - rflo)*(dflo - dval + 1) +
                             colr[rflo,dflo+1]*(rflo - rval + 1)*(dval - dflo) + 
                             colr[rflo+1,dflo+1]*(rval - rflo)*(dval - dflo))
    srcolors[tidx,bidx,2] = (colg[rflo,dflo]*(rflo - rval + 1)*(dflo - dval + 1) + 
                             colg[rflo+1,dflo]*(rval - rflo)*(dflo - dval + 1) +
                             colg[rflo,dflo+1]*(rflo - rval + 1)*(dval - dflo) + 
                             colg[rflo+1,dflo+1]*(rval - rflo)*(dval - dflo)) 
    srcolors[tidx,bidx,3] = (colb[rflo,dflo]*(rflo - rval + 1)*(dflo - dval + 1) + 
                             colb[rflo+1,dflo]*(rval - rflo)*(dflo - dval + 1) +
                             colb[rflo,dflo+1]*(rflo - rval + 1)*(dval - dflo) + 
                             colb[rflo+1,dflo+1]*(rval - rflo)*(dval - dflo))
    srcolors[tidx,bidx,4] = rrou
    srcolors[tidx,bidx,5] = drou
    
#%%****************************
#THIRD CELL - SAVE DATA TO FITS
#******************************       

greyscale_arr = (0.299*srcolors[:,:,1] + 0.587*srcolors[:,:,2] + 0.114*srcolors[:,:,3])
fits_arr = greyscale_arr.T #indexed by beta first then ejec_t

dustplotsave = os.path.join(imagedir, 'fitsdust')
if not os.path.exists(dustplotsave): os.makedirs(dustplotsave)
dustplotsave = os.path.join(dustplotsave, simin.split('\\')[-1][:-4])
dustplotfits = dustplotsave + '.fits'; dustplotvals = dustplotsave + '.txt'

if not os.path.exists(dustplotfits):
    hdu = fits.PrimaryHDU(fits_arr)    
    fitshdr = fits.Header()
    fitshdr['COMMENT'] = "Beta / Ejection time in file"
    hduhdr = fits.PrimaryHDU(header=fitshdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dustplotfits)

with open(dustplotvals, "w") as text_file:
    text_file.write('Ejection Times from ' + str(tvals[0]) + ' to ' +
    str(tvals[-1]) + ' days')
    text_file.write('\nBeta from ' + str(bvals[0]) + ' to ' +
    str(bvals[-1]))  
    text_file.write('\nEjection Time Values:\n ' +  str(tvals)[1:-1] +
                    '\nBeta Values:\n ' + str(bvals)[1:-1] )

#%%******************************
#SIXTH CELL - PLOT DATA ONTO IMAGE
#********************************

backgr_fill = (255,255,255,255)
featur_fill = (0,0,0,255)

logaxis = True
   
simtl = tvals[0]; simtu = tvals[-1]; betal = bvals[0]; betau = bvals[-1]

pixhi = 1000
pixwt = int(round(float(pixhi)/bno*tno))
border = 100
wscle = pixwt/(simtu - simtl)

if logaxis == True:
    hscle = pixhi/(np.log(betau) - np.log(betal))
else:
    hscle = pixhi/(betau - betal)        
    
dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
                             pixhi+int(3*border)),backgr_fill)
d = ImageDraw.Draw(dustimg)

imgmax = greyscale_arr.max()
if imgmax < 255: imgmax = 255

plotmethodlog = False
    
good_bool = srcolors[1:,1:,0] + srcolors[:-1,:-1,0] + srcolors[1:,:-1,0] + srcolors[:-1,1:,0]
good_locs = np.where(good_bool==4)

gaussed = gaussian_filter(greyscale_arr, sigma=20)
srdraw = greyscale_arr - gaussed

hih = 8
low = -8

fillco = np.clip(255.0/(hih-low)*(srdraw-low),0,255).astype(int)
for x in range(0, np.size(good_locs[0])):
    ta = good_locs[0][x]; ba = good_locs[1][x]
    logplotpixel(d,ta,ba,simres,dec,border,pixwt,pixhi,betal,simtl,
                  wscle,hscle,fillco[ta,ba],fillco[ta,ba],fillco[ta,ba])

a = d.polygon([(border,border),(border*2+pixwt,border), \
    (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
    outline = featur_fill)

[bminlocs,bmajlocs,tminlocs,tmajlocs,tmajticks,bmajticks,tminticks,bminticks] = bt_setaxisup(simtu,simtl,betau,betal,logaxis,border,pixhi,hscle,pixwt,wscle)  

majt = 20  #major tick length
mint = 10  #minor tick length
xaxis = pixhi + border*2
fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
fnt = ImageFont.truetype(fontloc, 20)
dtime = astropy.time.TimeDelta(1, format='jd')

for div in range(0, (np.size(bminlocs))): #beta axis minor ticks
    b = d.line([(border+mint,bminlocs[div]),(border,bminlocs[div])],\
    fill = featur_fill)
    
for div in range(0, (np.size(tminlocs))): #beta axis minor ticks
    b = d.line([(tminlocs[div],xaxis-mint),(tminlocs[div],xaxis)],\
    fill = featur_fill)

for div in range(0, (np.size(tmajlocs))): #simt axis major ticks
    b = d.line([(tmajlocs[div],xaxis-majt),(tmajlocs[div],xaxis)],\
    fill = featur_fill)
    ticktime = ctime - dtime*tmajticks[div]
    tick = ticktime.isot.replace('T','\n')[0:16]
    d.text((tmajlocs[div] - len(tick)*5,xaxis + 10), \
    tick, font=fnt, fill=featur_fill)

for div in range(0, (np.size(bmajlocs))): #beta axis major ticks
    b = d.line([(border+majt,bmajlocs[div]),(border,bmajlocs[div])],\
    fill = featur_fill)
    tick = "%1.1f" % bmajticks[div]
    d.text((border - len(tick)*5 - 40,bmajlocs[div] - 10 ), \
    tick, font=fnt, fill=featur_fill)
    
#axis labels
d.text((1.5*border + pixwt*0.5 - 150,pixhi + 2.7*border), \
"Date/Time of Ejection", font=fnt, fill=featur_fill)
d.text((0.25*border - 10,border-10), \
"Beta", font=fnt, fill=featur_fill)

#plot title
tfnt = ImageFont.truetype(fontloc, 30)
plttitle = (comdenom.upper() + ' ' + comname[:-1] + ' from ' + obsloc
+ '\n'+ ctime.isot[0:16].replace('T',' at '))
d.text((1.2*border,.25*border), \
plttitle, font=tfnt, fill=featur_fill)

dustimgfol = os.path.join(imagedir, 'dustplots')
if not os.path.exists(dustimgfol): os.makedirs(dustimgfol)
simparstr = simin.split('\\')[-1][:-4]
dustimgsave = os.path.join(dustimgfol,simparstr + '_gauss.png')
dustimg.save(dustimgsave,'png')
webbrowser.open(dustimgsave)