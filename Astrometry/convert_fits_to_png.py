import os
from astropy.io import fits
import numpy as np
from PIL import Image, ImageDraw


fitsdir = 'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1-MGN'
pngdir = 'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1-MGN\PNG_images'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".fts" in s]
fits_total = len(fits_list)

low = -0.7
hih = 1.25

for fits_id in range(0, fits_total):
        fitstemp = os.path.join(fitsdir, fits_list[fits_id])
        pngtemp = os.path.join(pngdir, (fits_list[fits_id].split(".")[0] +
        "_" + str(low) + "_" + str(hih) + ".png"))
        if not os.path.exists(pngtemp):
        
            hdulist = fits.open(fitstemp)
            data = (hdulist[0].data)
            img_size = np.shape(data)
            data = np.nan_to_num(data)
            
            pngimg = Image.new('RGBA', (img_size[0], img_size[1] ) , (0,0,0,255) )
            d = ImageDraw.Draw(pngimg)
                        
            filldat= np.clip(255.0/(hih-low)*(data-low),0,255).astype(int)            
            for x in xrange(0,img_size[0]):
                for y in xrange(0,img_size[1]):          
                    d.point((x,y), fill = (filldat[x,y],filldat[x,y],filldat[x,y]) )
            
    
            pngimg.save(pngtemp,'png')
        print fits_id