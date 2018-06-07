from astropy.io import fits
import os
import numpy as np
import sys
import matplotlib.pyplot as plt

#fits_folder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Dustplots\Stereo_A\HI-2-diff\NewFits'
#fits_out = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Dustplots\Stereo_A\HI-2-diff'

fits_folder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-2-diff'
fits_out = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-2-diff\star_filt'

fits_list = os.listdir(fits_folder)
fits_list = [s for s in fits_list if ".f" in s]

for fits_no in range(142,len(fits_list)):
    
    fits_string = fits_list[fits_no]
    fits_loc = os.path.join(fits_folder, fits_string)
    
    maskedimg = fits.open('C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\hi2A_mask.fts')
    imagemask  = maskedimg[0].data[::2,::2]

    data, header = fits.getdata(fits_loc, header=True)
    data = data*imagemask
    lim = 20000
    cdata = np.clip(data.flatten(),-lim,lim)
    cdata = cdata[cdata != 0]
    cdata = cdata[cdata != lim]
    cdata = cdata[cdata != -lim]
    #plt.hist(cdata,1000)
    
    standard = np.std(data)
    mean = np.mean(data)
    odata = np.copy(data)
    odata[np.where(data > (standard + mean))] = mean
    odata[np.where(data < -(standard + mean))] = mean
    
    fits_outfile = os.path.join(fits_out, fits_string)    
    fits.writeto(fits_outfile, odata, header , clobber=True)
    
    