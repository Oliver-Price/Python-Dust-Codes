import os
import numpy as np
from astropy.io import fits

backsave = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\ster_a_hi_1_background.fits'
sohosav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1'
subsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1\background subtracted'

img_list = sorted(os.listdir(sohosav))
img_list = [s for s in img_list if ".fts" in s]
img_total = len(img_list)

fitsimage = os.path.join(sohosav,img_list[0])
onedimg = fits.open(fitsimage)
values = onedimg[0].data
background = np.ones_like(values,dtype=float)
background.fill(99999999999)

for i in xrange(0,img_total):
    fitsimage = os.path.join(sohosav,img_list[i])
    onedimg = fits.open(fitsimage)
    values = onedimg[0].data
    values[np.isnan(values)]=99999999999
    background = np.minimum(background,values)
    
hdu = fits.PrimaryHDU(background)    
hdulist = fits.HDUList([hdu])
hdulist.writeto(backsave)

for i in xrange(0,img_total):
    fitsimage = os.path.join(sohosav,img_list[i])
    values, header_cur = fits.getdata(fitsimage, header=True)
    subimage = os.path.join(subsav,img_list[i])
    fits.writeto(subimage,(values-background),header_cur,clobber=True)
