from astropy.io import fits
import os
import numpy as np
import sys

fits_thresh_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue"
fits_wcs_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue"

fits_out_folder = os.path.join(fits_wcs_folder,'thresholded')
if not os.path.exists(fits_out_folder): os.makedirs(fits_out_folder)

fits_wcs_list = os.listdir(fits_wcs_folder)
fits_wcs_list = [s for s in fits_wcs_list if ".f" in s]

#npsave = os.path.join(os.path.dirname(fits_wcs_folder), os.path.basename(fits_wcs_folder) + '_thresholds.npy')
#nucvals = np.load(npsave)

#fitslistfile = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_B\good_astrometry.txt'
#threshfile = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_B\thresholdvalues.txt'

#with open(fitslistfile,'r') as dat:
#    fits_wcs_list = dat.readlines()
#fits_wcs_list = [f[:-1] for f in fits_wcs_list]
#    
#with open(threshfile,'r') as dat:
#    thresh_list = dat.readlines()
#thresh_list = [int(f) for f in thresh_list]
# 
coeff = 0.92 #0.97 stereo a hi-1, 0.92 soho

for fits_no in range(0,len(fits_wcs_list)):
    
    #fits_string = fits_wcs_list[fits_no][0:21]
    
    fits_wcs_loc = os.path.join(fits_wcs_folder, fits_wcs_list[fits_no])
    fits_thresh_loc = os.path.join(fits_thresh_folder, fits_wcs_list[fits_no])
    #fits_data_loc = os.path.join(fits_data_folder, fits_wcs_list[fits_no])
    
    data_wcs, header_wcs = fits.getdata(fits_wcs_loc, header=True)
    data_thresh, header_thresh = fits.getdata(fits_thresh_loc, header=True)
    #data_data, header_data = fits.getdata(fits_data_loc, header=True)
    
    #header_wcs['BITPIX'] = 32
    #header_wcs['NAXIS3'] = 1
    #data_wcs[np.where(data_thresh >= thresh_list[fits_no])] = 0
               
    '''
    if (np.size(nucvals) < len(fits_wcs_list)):
        data_wcs[np.where(data_wcs >= nucvals[sorted([0,fits_no,(np.size(nucvals)-1)])[1]]*coeff)] = 0
    else:
        data_wcs[np.where(data_wcs >= nucvals[fits_no]*coeff)] = 0
    '''
    data_wcs[np.where(data_wcs >= 13000)] = 0
    
    fits_outfile = os.path.join(fits_out_folder, fits_wcs_list[fits_no])    
    fits.writeto(fits_outfile, data_wcs, clobber=True)