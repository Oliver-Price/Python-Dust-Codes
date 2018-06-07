import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")

from conversion_routines import fixwraps

fitsdir = r'C:\PhD\Comet_data\Comet_NEAT_C2002V1\Gallery\Soho\C3_Clear_MGN'
rdcsav = r'C:\PhD\Comet_data\Comet_NEAT_C2002V1\Gallery\Soho\rdc.npy'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".fits" in s]
fits_total = len(fits_list)
radecorners = np.zeros((fits_total,8),dtype=float)

for fits_no in xrange(0, fits_total):

    fitsin = os.path.join(fitsdir, fits_list[fits_no])
    onedimg = fits.open(fitsin)
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
    print rafmax, rafmin
        
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
    
    print fits_no

np.save(rdcsav,radecorners)