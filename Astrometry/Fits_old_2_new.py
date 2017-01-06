from astropy.io import fits
import os
import numpy as np

fits_wcs_folder = r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\HI-1\calibrated pngs"
fits_data_folder = r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\HI-1\Badastrometry"
fits_out_folder = r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\HI-1\alignment_test"

fits_wcs_list = os.listdir(fits_wcs_folder)
fits_wcs_list = [s for s in fits_wcs_list if ".fits" in s]

for fits_no in xrange(0,len(fits_wcs_list)):
    
    fits_string = fits_wcs_list[fits_no][0:20] + 'B'
    
    fits_wcs_loc = os.path.join(fits_wcs_folder, fits_wcs_list[fits_no])
    fits_data_loc = os.path.join(fits_data_folder, fits_string + ".fts")
    
    data_wcs, header_wcs = fits.getdata(fits_wcs_loc, header=True)
    data_data, header_data = fits.getdata(fits_data_loc, header=True)
    
    #header_wcs['BITPIX'] = 32
    #header_wcs['NAXIS3'] = 1
    
    fits_outfile = os.path.join(fits_out_folder, fits_string + ".fits")    
    fits.writeto(fits_outfile, data_data.T, header_wcs, clobber=True)