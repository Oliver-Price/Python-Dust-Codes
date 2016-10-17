#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************

import string
import easygui
import os
import sys
import pickle
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")  

from orbitdata_loading_functions import orb_vector, orb_obs
from BT_plot_functions import beta2ypix, simt2xpix, plotpixel
from io_methods import get_obs_loc, get_stereo_instrument
from conversion_routines import fixwraps, round_to_base, RoundToSigFigs


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
if "Stereo" in obsloc: [sterinst, imagedir] = get_stereo_instrument(imagedir)
    
#import the orbit data
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,eq')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'eq,d10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

#choosing fits file to display and getting pathnames
pngdir = os.path.join(imagedir, 'cometplots')
pngin = easygui.fileopenbox(default = os.path.join(pngdir,'*'))
infile = os.path.basename(pngin)
filebase = infile[:string.find(infile,'.')]

#parameter savefile location
picklesavefile = os.path.join(pysav, 'imgsavs')
picklesavefile = os.path.join(picklesavefile, obsloc)
if "Stereo" in obsloc: picklesavefile = os.path.join(picklesavefile, sterinst)
picklesavefile = os.path.join(picklesavefile, filebase + '_plot_param.pickle')

greyscale = True

#check vital information exists
if not os.path.exists(picklesavefile):
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
if "Stereo" in obsloc: simsavbase = filebase[:filebase.find('A')+1]
else: simsavbase = filebase
simresdir = os.path.join(simresdir, obsloc)
simin = easygui.fileopenbox(default = os.path.join(simresdir, simsavbase + '*'))
simres = np.load(simin)

with open(simin[:-4] + '_parameters.pickle') as f:
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
    w = wcs.WCS(onedimg[0].header, key = 'A')
elif 'Soho' in obsloc:
    sys.exit("soho data not yet implemented")

#putting ra into sensible values
[ra_m, rafmin, rafmax] = fixwraps(ra, ramax, ramin)
    
#imagemask for HI-2
if "Stereo" in obsloc:
    if '2' in sterinst:
        maskedimg = fits.open('C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\hi2A_mask.fts')
        imagemask  = maskedimg[0].data[::2,::2]
        
    elif '1' in sterinst:
        direndindx = os.path.dirname(fitscoords).find("HI-1") + 4
        fits_dir = os.path.dirname(fitscoords)[:direndindx]
        fits_list = os.listdir(fits_dir)
        fitstr = [s for s in fits_list if filebase[:21] in s][0]
        orig_fits_file = os.path.join(fits_dir,fitstr)
        
        hdulist2 = fits.open(orig_fits_file)
        colours = (hdulist2[0].data)
        
        imagemask = np.ones_like(colr)
        imagemask[np.where(colours > 1.5e6)] = 0

else: imagemask = np.ones_like(colr)

#check if this has already been done
colmapsav = simin[:-4] + '_srcolors.npy'
locmapsav = simin[:-4] + '_pixelmapping.txt'

tno = np.size(tvals); bno = np.size(bvals)
srcolors = np.zeros((tno,bno,4),dtype=float)
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
    srcolors[tidx,bidx,0] = 1*imagemask[didx,ridx]
    srcolors[tidx,bidx,1] = colr[didx,ridx]
    srcolors[tidx,bidx,2] = colg[didx,ridx]
    srcolors[tidx,bidx,3] = colb[didx,ridx]

#%%****************************
#THIRD CELL - SAVE DATA TO FITS
#******************************       

greyscale_arr = (srcolors[:,:,1] + srcolors[:,:,2] + srcolors[:,:,3])/3
fits_arr = greyscale_arr.T #indexed by beta first then ejec_t

dustplotsave = os.path.join(imagedir, 'dustplots')
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
    simtl = tvals[0]; simtu = tvals[-1]; betal = bvals[0]; betau = bvals[-1]
    
    pixhi = 800
    pixwt = int(round(float(pixhi)/bno*tno))
    border = 100
    hscle = pixhi/(betau - betal)
    wscle = pixwt/(simtu - simtl)
        
    dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
                                 pixhi+int(3*border)),(0,0,0,255))
    d = ImageDraw.Draw(dustimg)
    
    imgmax = greyscale_arr.max()
    if imgmax < 255: imgmax = 255

    if 'Earth' in obsloc:
        plotmethodlog = False
    elif 'Stereo' in obsloc:
        if 'diff' or 'MGN' in sterinst: plotmethodlog = False
        else: plotmethodlog = True
    elif 'Soho' in obsloc:
        plotmethodlog = True
        sys.exit("soho data not yet implemented")

    if comdenom == 'c2011l4':
        low = 3000
        if obsloc == 'Stereo-B': hih = 20000
        elif obsloc == 'Stereo-A': hih = 70000
    elif comdenom == 'c2006p1':
        if "Stereo" in obsloc:
            low = 10000
            hih = 1500000
            if 'diff' in sterinst: 
                low = -1000
                hih = 1000
            elif 'MGN' in sterinst:  
                low = -0.7
                hih = 1.25
        elif obsloc == 'Earth':
            low = 0
            hih = 255
            
    good_bool = srcolors[1:,1:,0] + srcolors[:-1,:-1,0] + srcolors[1:,:-1,0] + srcolors[:-1,1:,0]
    good_locs = np.where(good_bool==4)

    if (plotmethodlog == True):
        
        grad = 1/(np.log10(hih) - np.log10(low))
        fillvals = np.clip(srcolors[:,:,1],1,9999999999)
        filcols = np.round(255*grad*(np.log10(fillvals) - np.log10(low)))
        filcols = np.clip(filcols,0,255).astype('int')
        
        for x in range(0, np.size(good_locs[0])):
            ta = good_locs[0][x]; ba = good_locs[1][x]
            plotpixel(d,ta,ba,simres,dec,border,pixwt,pixhi,betal,simtl,
                      wscle,hscle,filcols[ta,ba],filcols[ta,ba],filcols[ta,ba])
    else:
        fillco1 = np.clip(255.0/(hih-low)*(srcolors[:,:,1]-low),0,255).astype(int)
        fillco2 = np.clip(255.0/(hih-low)*(srcolors[:,:,2]-low),0,255).astype(int)
        fillco3 = np.clip(255.0/(hih-low)*(srcolors[:,:,3]-low),0,255).astype(int)
        for x in range(0, np.size(good_locs[0])):
            ta = good_locs[0][x]; ba = good_locs[1][x]
            plotpixel(d,ta,ba,simres,dec,border,pixwt,pixhi,betal,simtl,
                          wscle,hscle,fillco1[ta,ba],fillco2[ta,ba],fillco3[ta,ba])
            
#%%**********************
#SEVENTH CELL - DRAW AXIS
#************************
    
    a = d.polygon([(border,border),(border*2+pixwt,border), \
        (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
        outline = (255,255,255,128))
    
    bl1sf = RoundToSigFigs(betal,1); bu1sf = RoundToSigFigs(betau,1)
    tl1sf = round(simtl); tu1sf = round(simtu)
    
    tdivmajors = np.array([1.,2.,3.])
    tdivnos = (1/tdivmajors) * (tu1sf - tl1sf)
    tnodivs = 6
    tdividx = np.where((tdivnos <= tnodivs)==True)[0][0]
    tmajticks = np.arange(tu1sf, tl1sf-0.0001, -tdivmajors[tdividx])
    
    tdivminors = np.array([0.5,1,1])
    tminticks = np.arange(tu1sf, tl1sf-0.0001, -tdivminors[tdividx])
    tminticks = np.setdiff1d(tminticks,tmajticks)
            
    bdivmajors = np.array([0.1,0.2,0.5,1,2])
    bdivnos = (1/bdivmajors) * (bu1sf - bl1sf)
    bnodivs = 10
    bdividx = np.where((bdivnos <= bnodivs)==True)[0][0]
    bu2majdv = round_to_base(betau, bdivmajors[bdividx])
    bl2majdv = round_to_base(betal, bdivmajors[bdividx])   
    bmajticks = np.arange(bu2majdv, bl2majdv-0.0001, -bdivmajors[bdividx])
    
    bdivminors = np.array([0.02,0.05,0.1,0.2,0.5])
    bu2mindv = round_to_base(betau, bdivminors[bdividx])
    bl2mindv = round_to_base(betal, bdivminors[bdividx])   
    bminticks = np.arange(bu2mindv, bl2mindv-0.0001, -bdivminors[bdividx])
    bminticks = np.setdiff1d(bminticks,bmajticks)
    
    bminlocs = beta2ypix(bminticks, border, pixhi, betal, hscle)
    tminlocs = simt2xpix(tminticks, border, pixwt, simtl, wscle)
    bmajlocs = beta2ypix(bmajticks, border, pixhi, betal, hscle)
    tmajlocs = simt2xpix(tmajticks, border, pixwt, simtl, wscle)
    
    majt = 20  #major tick length
    mint = 10  #minor tick length
    xaxis = pixhi + border*2
    fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
    fnt = ImageFont.truetype(fontloc, 20)
    dtime = astropy.time.TimeDelta(1, format='jd')
    
    for div in xrange(0, (np.size(bminlocs))): #beta axis minor ticks
        b = d.line([(border+mint,bminlocs[div]),(border,bminlocs[div])],\
        fill = (255,255,255,128))
        
    for div in xrange(0, (np.size(tminlocs))): #beta axis minor ticks
        b = d.line([(tminlocs[div],xaxis-mint),(tminlocs[div],xaxis)],\
        fill = (255,255,255,128))
    
    for div in xrange(0, (np.size(tmajlocs))): #simt axis major ticks
        b = d.line([(tmajlocs[div],xaxis-majt),(tmajlocs[div],xaxis)],\
        fill = (255,255,255,128))
        ticktime = ctime - dtime*tmajticks[div]
        tick = string.replace(ticktime.isot,'T','\n')[0:16]
        d.text((tmajlocs[div] - len(tick)*5,xaxis + 10), \
        tick, font=fnt, fill=(255,255,255,128))
    
    for div in xrange(0, (np.size(bmajlocs))): #beta axis major ticks
        b = d.line([(border+majt,bmajlocs[div]),(border,bmajlocs[div])],\
        fill = (255,255,255,128))
        tick = str(bmajticks[div])
        d.text((border - len(tick)*5 - 40,bmajlocs[div] - 10 ), \
        tick, font=fnt, fill=(255,255,255,128))
        
    #axis labels
    d.text((1.5*border + pixwt*0.5 - 150,pixhi + 2.7*border), \
    "Date/Time of Ejection", font=fnt, fill=(255,255,255,128))
    d.text((0.25*border - 10,border-10), \
    "Beta", font=fnt, fill=(255,255,255,128))
    
    #plot title
    tfnt = ImageFont.truetype(fontloc, 30)
    plttitle = (comdenom.upper() + ' ' + comname[:-1] + ' from ' + obsloc
    + '\n'+ string.replace(ctime.isot[0:16],'T',' at '))
    d.text((1.2*border,.25*border), \
    plttitle, font=tfnt, fill=(255,255,255,128))
    
    dustimg.show()
    dustimgsav = simin[:-4] + '.png'  
    dustimg.save(dustimgsav,'png')