import os
import numpy as np
from astropy.io import fits

backsave = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\soho_clear_background.fits'
sohosav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\processing\Renamed Data'
#subsav = os.path.join(sohosav,'background_subtracted')
subsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\processing\background subtracted'
if not os.path.exists(subsav): os.makedirs(subsav)


img_list = sorted(os.listdir(sohosav))
img_list = [s for s in img_list if ".fits" in s]
img_total = len(img_list)

fitsimage = os.path.join(sohosav,img_list[0])
onedimg = fits.open(fitsimage)
values = onedimg[0].data
background = np.ones_like(values,dtype=float)
background.fill(99999999999)

for i in range(0,img_total):
    fitsimage = os.path.join(sohosav,img_list[i])
    onedimg = fits.open(fitsimage)
    values = onedimg[0].data
    values[np.isnan(values)]=values.max()
    background = np.minimum(background,values)
    
hdu = fits.PrimaryHDU(background)    
hdulist = fits.HDUList([hdu])
hdulist.writeto(backsave)

for i in range(0,img_total):
    fitsimage = os.path.join(sohosav,img_list[i])
    values, header_cur = fits.getdata(fitsimage, header=True)
    subimage = os.path.join(subsav,img_list[i])
    fits.writeto(subimage,(values-background),clobber=True)
