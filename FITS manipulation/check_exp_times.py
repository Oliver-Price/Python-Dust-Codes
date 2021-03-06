# -*- coding: utf-8 -*-
#Created on Tue Feb 14 11:58:37 2017

from astropy.io import fits
import os
import numpy as np

fits_in_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Raw"
fits_out_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Uncalibrated"

expsave = 'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\exptimes.npy'

fits_in_list = os.listdir(fits_in_folder)
fits_in_list = [s for s in fits_in_list if ".fts" in s]
fits_out_list = []
fits_no_list = []

exptimes = np.zeros(shape=(196))

for fits_no in xrange(0,len(fits_in_list)):
    
    fits_in_loc = os.path.join(fits_in_folder, fits_in_list[fits_no])
    data_raw, header_raw = fits.getdata(fits_in_loc, header=True)
    
    fits_new_name = header_raw['DATE-OBS'].replace('/','') + '_' + header_raw['TIME-OBS'][0:8].replace(':','') + '_Clear.fits'
    
    fits_out_list.append(fits_new_name)     
    exptimes[fits_no] = np.round(header_raw['EXPTIME'],1)
    fits_no_list.append(fits_in_list[fits_no][4:8])
    
    print np.round(header_raw['EXPTIME'],2), fits_new_name, fits_in_list[fits_no][4:8]
