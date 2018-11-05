import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
import easygui

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from io_methods import get_obs_loc, get_stereo_instrument, get_soho_instrument
from orbitdata_loading_functions import orb_vector, orb_obs
from conversion_routines import fixwraps

rdcsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\rdc-HI-1.npy'
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

fitsdir=imagedir

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)
radecorners = np.zeros((fits_total,8),dtype=float)

for fits_no in range(0, fits_total):

    fitsin = os.path.join(fitsdir, fits_list[fits_no])
    onedimg = fits.open(fitsin)
    try:
        w = wcs.WCS(onedimg[0].header, key = 'A')
    except:
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
    
    [ra_m, rafmin, rafmax, bool_val] = fixwraps(ra, ramax, ramin)
        
    try:
        if rafmax > trafmax:
            trafmax = rafmax
    except: trafmax = rafmax
    
    try:
        if rafmin < trafmin:
            trafmin = rafmin
    except: trafmin = rafmin  
    
    try:
        if decmax > tdecmax:
            tdecmax = decmax
    except: tdecmax = decmax
    
    try:
        if decmin < tdecmin:
            tdecmin = decmin
    except: tdecmin = decmin
    
    radecorners[fits_no,0] = ra[0,0]; radecorners[fits_no,1] = ra[-1,0]
    radecorners[fits_no,2] = ra[0,-1]; radecorners[fits_no,3] = ra[-1,-1]
    radecorners[fits_no,4] = dec[0,0]; radecorners[fits_no,5] = dec[-1,0]
    radecorners[fits_no,6] = dec[0,-1]; radecorners[fits_no,7] = dec[-1,-1]
    
    print(fits_no)

np.save(rdcsav,radecorners)