# -*- coding: utf-8 -*-
#Created on Tue Feb 14 11:58:37 2017

from astropy.io import fits
import os
import numpy as np

fits_in_folder = r"C:\Users\op2\Downloads\86798720"
fits_out_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Orange"

#expsave = 'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\orangeexptimes.npy'

fits_in_list = os.listdir(fits_in_folder)
fits_in_list = [s for s in fits_in_list if ".f" in s]

#exptimes = np.zeros(shape=(129))

for fits_no in range(0,len(fits_in_list)):
    
    fits_in_loc = os.path.join(fits_in_folder, fits_in_list[fits_no])
    data_raw, header_raw = fits.getdata(fits_in_loc, header=True)
    
    fits_new_name = (header_raw['DATE-OBS'].replace('/','') + '_' + header_raw['TIME-OBS'][0:8].replace(':','') + '_Orange.fits')
    #fits_new_name = header_raw['DATE-OBS'][:19].replace('-','').replace(':','').replace('T','_') + '_Bl.fits'
    
    fits_outfile = os.path.join(fits_out_folder, fits_new_name)     
    #exptimes[fits_no] = np.round(header_raw['EXPTIME'],1)

    header = header_raw[:77]
    fits.writeto(fits_outfile, data_raw, header, clobber=True)