import os
from astropy.io import fits
import numpy as np
from PIL import Image, ImageDraw


fitsdir = 'C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B'
pngdir = 'C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\pngs'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".fts" in s]
fits_total = len(fits_list)

low = 10000
hih = 1500000

for fits_id in range(0, fits_total):
        fitstemp = os.path.join(fitsdir, fits_list[fits_id])
        pngtemp = os.path.join(pngdir, (fits_list[fits_id].split(".")[0] + ".png"))
        if not os.path.exists(pngtemp):
        
            hdulist = fits.open(fitstemp)
            data = (hdulist[0].data)
            img_size = np.shape(data)
            data = np.nan_to_num(data)
            
            pngimg = Image.new('RGBA', (img_size[0], img_size[1] ) , (0,0,0,255) )
            d = ImageDraw.Draw(pngimg)
            
            for x in xrange(0,img_size[0]):
                for y in xrange(0,img_size[1]):
                    
                    fillval = sorted([1, data[x,y], 9999999999])[1]
                    fillco = int(round(255*(np.log10(fillval) - np.log10(low))*
    								1 / (np.log10(hih) - np.log10(low))))
                    fillco = sorted([0, fillco, 255])[1]
            
                    d.point( (x,y) , fill = (fillco,fillco,fillco) )
            
    
            pngimg.save(pngtemp,'png')
        print fits_id