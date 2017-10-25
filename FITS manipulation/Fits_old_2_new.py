from astropy.io import fits
import os
import numpy as np
import sys

#fits_wcs_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\Lasco Processing\Exptime_Sorted\8.1\background_calibrated"
#fits_data_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\Lasco Processing\Uncalibrated"
#fits_out_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\Soho\C3_Clear"

fits_wcs_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\processing\Astrometry_headers_bckgr"
fits_out_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear"
fits_data_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\processing\Renamed Data"

if not os.path.exists(fits_out_folder): os.makedirs(fits_out_folder)

fits_wcs_list = os.listdir(fits_wcs_folder)
fits_wcs_list = [s for s in fits_wcs_list if ".f" in s]

for fits_no in range(0,4):#len(fits_wcs_list)):
    
    fits_string = fits_wcs_list[fits_no][0:21]
    
    fits_wcs_loc = os.path.join(fits_wcs_folder, fits_wcs_list[fits_no])
    fits_data_loc = os.path.join(fits_data_folder, fits_wcs_list[fits_no])

    data_wcs, header_wcs = fits.getdata(fits_wcs_loc, header=True)
    data_data, header_data = fits.getdata(fits_data_loc, header=True)
    
    #header_wcs['BITPIX'] = 32
    #header_wcs['NAXIS3'] = 1
    
    fits_outfile = os.path.join(fits_out_folder, fits_wcs_list[fits_no])    
    fits.writeto(fits_outfile, data_data, header_wcs , clobber=True)