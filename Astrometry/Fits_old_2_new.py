from astropy.io import fits
import os
import numpy as np

fits_wcs_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\processing\Astrometry_headers_bckgr"
fits_data_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear_MGN\oldastro"
fits_out_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear_MGN"

fits_wcs_list = os.listdir(fits_wcs_folder)

for fits_no in xrange(0,len(fits_wcs_list)):
    
    fits_string = fits_wcs_list[fits_no][0:21]
    
    fits_wcs_loc = os.path.join(fits_wcs_folder, fits_wcs_list[fits_no])
    fits_data_loc = os.path.join(fits_data_folder, fits_string + "_mgn.fits")
    
    data_wcs, header_wcs = fits.getdata(fits_wcs_loc, header=True)
    data_data, header_data = fits.getdata(fits_data_loc, header=True)
    
    #header_wcs['BITPIX'] = 32
    #header_wcs['NAXIS3'] = 1
    
    fits_outfile = os.path.join(fits_out_folder, fits_string + "_mgn.fits")    
    fits.writeto(fits_outfile, data_data, header_wcs, clobber=True)