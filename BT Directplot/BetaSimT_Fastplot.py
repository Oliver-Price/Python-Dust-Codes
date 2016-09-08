#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************

import string
import easygui
import os
import sys
import pickle
import numpy as np
import time
import astropy
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")  

from orbitdata_loading_functions import orb_vector, orb_obs
from BT_plot_functions import beta2ypix, linsimt2xpix, logsimt2xpix, radec_slim, \
greyscale_remap


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
bool_locs = np.array([(('EARTH' in obslocstr) or ('Earth' in obslocstr)),
                 (('STEREO-A' in obslocstr) or ('Stereo-A' in obslocstr)
                 or ('STEREO_A' in obslocstr) or ('Stereo_A' in obslocstr)),
                 (('STEREO-B' in obslocstr) or ('Stereo-B' in obslocstr)
                 or ('STEREO_B' in obslocstr) or ('Stereo_B' in obslocstr)),
                 (('SOHO' in obslocstr) or ('Soho' in obslocstr))])
name_locs = np.array(['Earth', 'Stereo_A', 'Stereo_B', 'Soho'])
case_locs = np.size(np.nonzero(bool_locs))
if case_locs > 1:
    obsmsg = "Please select observer location"
    obschoices = name_locs[bool_locs].tolist()
    obsloc = easygui.buttonbox(obsmsg, choices=obschoices)
    imagedir = os.path.join(imagedir, obsloc)
elif case_locs == 1:
    obsloc = name_locs[bool_locs][0]
else: sys.exit("No Good Observer Location")
    
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

#parameter savefile location
picklesavefile = os.path.join(pysav, filebase + '_dustplot')
picklexists = os.path.exists(picklesavefile)

greyscale = True

#check vital information exists
if (picklexists == False):
    sys.exit("Image not calculated")

#import important information
with open(picklesavefile) as f:
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
simin = easygui.fileopenbox(default = os.path.join(simresdir, filebase + '*'))
simres = np.load(simin)

with open(simin[:-4] + '_parameters') as f:
    sparameters = pickle.load(f)
    tmax = sparameters[0]
    tno = sparameters[2]
    bno = sparameters[3]
    tspace = sparameters[4]

#%%*********************************************************
#SECOND CELL - IDENTIFYING PIXEL VALUES TO USE FOR BETA/SIMT
#***********************************************************

onedimg = fits.open(fitscoords)
if 'Earth' in obsloc:
    w = wcs.WCS(onedimg[0].header)
elif 'Stereo' in obsloc:
    w = wcs.WCS(onedimg[0].header, key = 'A')
elif 'Soho' in obsloc:
    sys.exit("soho data not yet implemented")

#putting ra into sensible values
if ramax < 0:
    ra_m = np.copy(ra) + 360
    rafmin = np.amin(ra_m)
    rafmax = np.amax(ra_m)
else:
    ra_m = ra
    rafmin = ramin
    rafmax = ramax

#check if this has already been done
colmapsav = simin[:-4] + '_srcolors.npy'
locmapsav = simin[:-4] + '_pixelmapping.txt'
colmapsavexists = os.path.exists(colmapsav)
a = time.clock()

if (colmapsavexists == True): #load if it has
    srcolors = np.load(colmapsav)
    
elif (colmapsavexists == False): #do if it hasnt
    srcolors = np.zeros((tno,bno,4),dtype=int)
    simres_rounded = np.round(simres[:,:,14:17]).astype(int)
    simres_rounded = np.clip(simres_rounded, 0, 1023)    
    
    non_zeros_0 = simres_rounded[:,:,0].nonzero()[0]
    non_zeros_1 = simres_rounded[:,:,0].nonzero()[1]
    no_points = non_zeros_0.size
    
    for x in xrange(0,no_points):
        tidx = non_zeros_0[x]
        bidx = non_zeros_1[x]
        ridx = simres_rounded[tidx,bidx,1]
        didx = simres_rounded[tidx,bidx,2]
        srcolors[tidx,bidx,0] = 1
        srcolors[tidx,bidx,1] = colr[didx,ridx]
        srcolors[tidx,bidx,2] = colg[didx,ridx]
        srcolors[tidx,bidx,3] = colb[didx,ridx]

b = time.clock()
#%%****************************
#THIRD CELL - SAVE DATA TO FITS
#******************************       

greyscale_arr = (srcolors[:,:,1] + srcolors[:,:,2] + srcolors[:,:,3])/3
fits_arr = greyscale_arr.T #indexed by beta first then ejec_t

dustplotsave = os.path.join(imagedir, 'dustplots')
if not os.path.exists(dustplotsave):
    os.makedirs(dustplotsave)
dustplotsave = os.path.join(dustplotsave, simin.split('\\')[-1][:-4])
dustplotfits = dustplotsave + '.fits'
dustplotvals = dustplotsave + '.txt'

if not os.path.exists(dustplotfits):
    hdu = fits.PrimaryHDU(fits_arr)    
    fitshdr = fits.Header()
    fitshdr['COMMENT'] = "Beta / Ejection time in file"
    hduhdr = fits.PrimaryHDU(header=fitshdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dustplotfits)

with open(dustplotvals, "w") as text_file:
    text_file.write('Ejection Times from ' + str(simres[0,0,0]) + ' to ' +
    str(simres[-1,0,0]) + ' days, ')
    if tspace == 'Linear':
        text_file.write('with a linear spacing of ' + 
                        str(np.round(24*60*(simres[1,0,0] - simres[0,0,0]),1))
                        + ' minutes')  
    text_file.write('\nBeta from ' + str(simres[0,0,1]) + ' to ' +
    str(simres[-1,0,0]))  
    text_file.write('\nEjection Time Values:\n ' +  str(simres[:,0,0])[1:-1] +
                    '\nBeta Values:\n ' + str(simres[0,:,1])[1:-1] )

#%%**************************************
#FOURTH CELL - SAVE TIME FILTERED TO FITS
#****************************************      
filtermsg = "Create ejection time filtered fits?"
reply = easygui.ynbox(msg=filtermsg)

if reply == True:
    
    hmsg = "Choose time filter shift in hours"
    hreply = easygui.enterbox(hmsg)
    while 1:
        if hreply == None: break #exit if cancel pressed                 
        errmsg = ""
        try:
            tshift = float(hreply)
        except ValueError:
                errmsg += ("Time must be a number")
        if errmsg == "": break
        hreply = easygui.enterbox(errmsg)
    tshift = int(float(hreply)/(simres[1,0,0] - simres[0,0,0])/24)
    
    dustplotmodifits = (dustplotsave + '_tfilter_' +
                        string.replace(str(float(hreply)),'.','\'') + '.fits') 
    
    if not os.path.exists(dustplotmodifits):
        ref_vals = np.copy(simres[:,:,15].T)
        
        fits_t_minus = np.copy(fits_arr[:,tshift:])
        fits_t_plus = np.copy(fits_arr[:,:-tshift])
        ref_t_minus = np.copy(ref_vals[:,tshift:])
        ref_t_plus = np.copy(ref_vals[:,:-tshift])        
        
        fits_t_modified = np.copy(2*fits_arr)
        ref_t_modified = np.copy(ref_vals)
        fits_t_modified[:,:-tshift] = fits_t_modified[:,:-tshift] - fits_t_minus
        fits_t_modified[:,tshift:] = fits_t_modified[:,tshift:] - fits_t_plus
        ref_t_modified[:,:-tshift] = ref_t_modified[:,:-tshift] + ref_t_minus
        ref_t_modified[:,tshift:] = ref_t_modified[:,tshift:] + ref_t_plus
        
        fits_t_modified[:,0:tshift+1] = 0
        fits_t_modified[:,-tshift-1:] = 0 
        fits_t_modified[np.where(ref_t_modified!=3)] = 0
        fits_t_modified[np.where(fits_t_modified<0)] = 0
        #fits_t_modified[np.where(fits_t_modified>255)] = 0
        
        #fits_t_modified[np.where(fits_t_modified>50)] = 0 #EXPERIMENTAL
        
        hdu = fits.PrimaryHDU(fits_t_modified)    
        fitshdr = fits.Header()
        fitshdr['COMMENT'] = "Beta / Ejection time in file"
        fitshdr['COMMENT'] = "Modified with a Larson-Sekanina(esque) filter"
        hduhdr = fits.PrimaryHDU(header=fitshdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dustplotmodifits)
        
#%%*************************************
#FIFTH CELL - SAVE BETA FILTERED TO FITS
#***************************************
filtermsg = "Create beta filtered fits?"
reply = easygui.ynbox(msg=filtermsg)

if reply == True:
    
    hmsg = "Choose beta filter shift (cell no)"
    hreply = easygui.enterbox(hmsg)
    while 1:
        if hreply == None: break #exit if cancel pressed                 
        errmsg = ""
        try:
            bshift = int(hreply)
        except ValueError:
                errmsg += ("Beta must be a number")
        if errmsg == "": break
        hreply = easygui.enterbox(errmsg)
    bshift = int(hreply)
    
    fits_b_minus = fits_arr[bshift:,:]
    fits_b_plus = fits_arr[:-bshift,:]
    
    fits_b_modified = 2*fits_arr
    fits_b_modified[:-bshift,:] = fits_b_modified[:-bshift,:] - fits_b_minus
    fits_b_modified[bshift:,:] = fits_b_modified[bshift:,:] - fits_b_plus
    
    fits_b_modified[np.where(fits_b_modified<0)] = 0
    #fits_b_modified[np.where(fits_b_modified>255)] = 255
    
    dustplotmodifits = (dustplotsave + '_bfilter_' +
                        string.replace(str(float(hreply)),'.','\'') + '.fits')
    
    if not os.path.exists(dustplotmodifits):
        hdu = fits.PrimaryHDU(fits_b_modified)    
        fitshdr = fits.Header()
        fitshdr['COMMENT'] = "Beta / Ejection time in file"
        fitshdr['COMMENT'] = "Modified with a Larson-Sekanina(esque) filter"
        hduhdr = fits.PrimaryHDU(header=fitshdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dustplotmodifits)
                  
#%%******************************
#SIXTH CELL - PLOT DATA ONTO IMAGE
#*********************************          
 
#choose img type
colormsg = "Preview dustplot output?"
reply = easygui.ynbox(msg=colormsg)

if reply == True:    
    if tspace == 'Logarithmic':
        sys.exit('Only a linear timescale can be displayed')     
        
    simtl = simres[0,0,0]; simtu = simres[tno-1,0,0]
    betal = simres[0,0,1]; betau = simres[0,bno-1,1]
    
    t2sfu = float('%.2g' % simtu)
    t2sfl = float('%.2g' % simtl)
    b1sfu = float('%.1g' % betau)
    b1sfl = float('%.1g' % betal)
    
    pixhi = 1200
    pixwt = int(round(float(pixhi)/bno*tno))
    border = 100
    hscle = pixhi/(np.log10(betau) - np.log10(betal))
    wscle = pixwt/(simtu - simtl)
        
    dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
                                 pixhi+int(3*border)),(0,0,0,255))
    d = ImageDraw.Draw(dustimg)
    nmmax = np.log(np.max(srcolors[:,:,3])+1)
    
    imgmax = greyscale_arr.max()
    if imgmax < 255: imgmax = 255
        
    fudgefactor = 0.8
    if comdenom == 'c2011l4':
        low = 3000
        hih = 20000
    elif comdenom == 'c2006p1':
        low = 10000
        hih = 1500000  
    #newmap = greyscale_remap(200,50,mode = 'Linear')
        
    greyscale_disp = True
    if (greyscale_disp == True):  
        for ta in xrange(0, tno-1):
            for ba in xrange(0, bno-1):
                fillval = sorted([1, greyscale_arr[ta,ba], 9999999999])[1]
                fillco = int(round(255*(np.log10(fillval) - np.log10(low))*
                					1 / (np.log10(hih) - np.log10(low))))
                fillco = int(round(greyscale_arr[ta,ba]*255/imgmax/fudgefactor))
                fillco = sorted([0, fillco, 255])[1]
#               fillco = greyscale_arr[ta,ba]
                b1 = beta2ypix(simres[ta,ba,1], border, pixhi, b1sfl, hscle)
                t1 = linsimt2xpix(simres[ta,ba,0], border, t2sfl, wscle)
                b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, b1sfl, hscle)
                t2 = linsimt2xpix(simres[ta,ba+1,0], border, t2sfl, wscle)
                b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, b1sfl, hscle)
                t3 = linsimt2xpix(simres[ta+1,ba+1,0], border, t2sfl, wscle)
                b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, b1sfl, hscle)
                t4 = linsimt2xpix(simres[ta+1,ba,0], border, t2sfl, wscle)
                a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
                ,fill=(fillco,fillco,fillco,255))
    else:
        for ta in xrange(0, tno-1):
            for ba in xrange(0, bno-1):
                b1 = beta2ypix(simres[ta,ba,1], border, pixhi, b1sfl, hscle)
                t1 = linsimt2xpix(simres[ta,ba,0], border, t2sfl, wscle)
                b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, b1sfl, hscle)
                t2 = linsimt2xpix(simres[ta,ba+1,0], border, t2sfl, wscle)
                b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, b1sfl, hscle)
                t3 = linsimt2xpix(simres[ta+1,ba+1,0], border, t2sfl, wscle)
                b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, b1sfl, hscle)
                t4 = linsimt2xpix(simres[ta+1,ba,0], border, t2sfl, wscle)
                a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
                ,fill=(srcolors[ta,ba,0],srcolors[ta,ba,1],srcolors[ta,ba,2],255))
                
#%%**********************
#SEVENTH CELL - DRAW AXIS
#************************
    
    a = d.polygon([(border,border),(border*2+pixwt,border), \
        (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
        outline = (255,255,255,128))
    
    decades = np.logspace(-4,4,9)
    
    bdecl = np.searchsorted(decades,betal, side = 'right')-1
    bdecu = np.searchsorted(decades,betau)
    
    bminticks = np.linspace(decades[bdecl],decades[bdecl+1],10)
    for bdec in xrange(bdecl+1, bdecu):
        bminticks = np.concatenate((bminticks,
                          np.linspace(decades[bdec],decades[bdec+1],10)[1:10]))
    bminticks = bminticks[np.searchsorted(bminticks,b1sfl):
                              np.searchsorted(bminticks,b1sfu)+1]
    bmajticks = np.intersect1d(bminticks,decades)
    bmajticks = np.union1d(bmajticks,np.array([b1sfl]))
    bmajticks = np.union1d(bmajticks,np.array([b1sfu]))
    bminticlocs = beta2ypix(bminticks, border, pixhi, b1sfl, hscle)
    bmajticlocs = beta2ypix(bmajticks, border, pixhi, b1sfl, hscle)
    
    lindivmajors = np.array([0.1,0.2,0.5,1,2,5,10,20,50,100,200])
    lindivrecips = np.array([10,5,2,1,0.5,0.2,0.1,0.05,0.02,0.01,0.005])
    
    tdivnos = lindivrecips * (t2sfu - t2sfl)
    nodivs = 6
    tdividx = (np.abs(tdivnos-nodivs)).argmin()
    tlodi = np.floor(t2sfl*lindivrecips[tdividx])*lindivmajors[tdividx]
    thidi = np.ceil(t2sfu*lindivrecips[tdividx])*lindivmajors[tdividx]
    tmajticks = np.arange(tlodi, thidi+1e-10, lindivmajors[tdividx])
        
    #tminticlocs = linsimt2xpix(tminticks, border, t1sfl, wscle)
    tmajticlocs = linsimt2xpix(tmajticks, border, t2sfl, wscle)
    
    majt = 20  #major tick length
    mint = 10  #minor tick length
    xaxis = pixhi + border*2
    fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
    fnt = ImageFont.truetype(fontloc, 20)
    dtime = astropy.time.TimeDelta(1, format='jd')
    
    for div in xrange(0, (np.size(bminticlocs))): #beta axis minor ticks
        b = d.line([(border+mint,bminticlocs[div]),(border,bminticlocs[div])],\
        fill = (255,255,255,128))
    
    for div in xrange(0, (np.size(tmajticlocs))): #simt axis major ticks
        b = d.line([(tmajticlocs[div],xaxis-majt),(tmajticlocs[div],xaxis)],\
        fill = (255,255,255,128))
        ticktime = ctime - dtime*tmajticks[div]
        tick = string.replace(ticktime.isot,'T','\n')[0:16]
        d.text((tmajticlocs[div] - len(tick)*5,xaxis + 10), \
        tick, font=fnt, fill=(255,255,255,128))
    
    for div in xrange(0, (np.size(bmajticlocs))): #beta axis major ticks
        b = d.line([(border+majt,bmajticlocs[div]),(border,bmajticlocs[div])],\
        fill = (255,255,255,128))
        tick = str(bmajticks[div])
        d.text((border - len(tick)*5 - 40,bmajticlocs[div] - 10 ), \
        tick, font=fnt, fill=(255,255,255,128))
        
    #axis labels
    d.text((1.5*border + pixwt*0.5 - 150,pixhi + 2.7*border), \
    "Date/Time of Ejection", font=fnt, fill=(255,255,255,128))
    d.text((0.25*border - 10,0.75*border - 20), \
    "Beta", font=fnt, fill=(255,255,255,128))
    
    #plot title
    plttitle = (comdenom.upper() + ' ' + comname[:-1] + ' from ' + 
    obsloc + ' from date: ' + string.replace(ctime.isot[0:16],'T',' at '))
    tfnt = ImageFont.truetype(fontloc, 30)
    d.text((1.5*border + pixwt*0.5 - len(plttitle)*5 - 200,.35*border), \
    plttitle, font=tfnt, fill=(255,255,255,128))
    
    dustimg.show()
    dustimgsav = simin[:-4] + '.png'  
    dustimg.save(dustimgsav,'png')