# -*- coding: utf-8 -*-
#Created on Tue Feb 14 14:20:15 2017

from astropy.io import fits
import os
import numpy as np

fits_in_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Raw"
fits_out_folder = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Uncalibrated"

expsave = 'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\exptimes.npy'

fits_in_list = os.listdir(fits_in_folder)
fits_in_list = [s for s in fits_in_list if ".fts" in s]
fits_name_list = []

for fits_no in xrange(0,len(fits_in_list)):
    
    fits_in_loc = os.path.join(fits_in_folder, fits_in_list[fits_no])
    data_raw, header_raw = fits.getdata(fits_in_loc, header=True)
    fits_new_name = header_raw['DATE-OBS'].replace('/','') + '_' + header_raw['TIME-OBS'][0:8].replace(':','') + '_Clear.fits'
    fits_name_list.append(fits_new_name)
    
exptimes = np.load(expsave)
fits_list = np.array(fits_name_list)

fits_4p1 = fits_list[np.where(exptimes == 4.1)]
fits_8p1 = fits_list[np.where(exptimes == 8.1)]
fits_19p1 = fits_list[np.where(exptimes == 19.1)]

fits_leftovers = np.setdiff1d(fits_list,fits_4p1)
fits_leftovers = np.setdiff1d(fits_leftovers,fits_8p1)
fits_leftovers = np.setdiff1d(fits_leftovers,fits_19p1)

fits_4p1_save = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\4.1"
fits_8p1_save = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\8.1"
fits_19p1_save = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\19.1"
fits_others_save = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\Others"
fits_others_bgkr = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\Others\background_subtracted"

for fits_no in xrange(0,len(fits_4p1)):
    fits_in_loc = os.path.join(fits_out_folder, fits_4p1[fits_no])
    fits_out_loc = os.path.join(fits_4p1_save, fits_4p1[fits_no])
    data_raw, header_raw = fits.getdata(fits_in_loc, header=True)
    fits.writeto(fits_out_loc, data_raw, clobber=True)
    
for fits_no in xrange(0,len(fits_8p1)):
    fits_in_loc = os.path.join(fits_out_folder, fits_8p1[fits_no])
    fits_out_loc = os.path.join(fits_8p1_save, fits_8p1[fits_no])
    data_raw, header_raw = fits.getdata(fits_in_loc, header=True)
    fits.writeto(fits_out_loc, data_raw, clobber=True)
    
for fits_no in xrange(0,len(fits_19p1)):
    fits_in_loc = os.path.join(fits_out_folder, fits_19p1[fits_no])
    fits_out_loc = os.path.join(fits_19p1_save, fits_19p1[fits_no])
    data_raw, header_raw = fits.getdata(fits_in_loc, header=True)
    fits.writeto(fits_out_loc, data_raw, clobber=True)

backfolder = r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\background_guesses'
back_8p1 = fits.getdata(r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\soho_8.1_background.fits')
back_4p1 = fits.getdata(r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\soho_4.1_background.fits')
back_19p1 = fits.getdata(r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\soho_19.1_background.fits')

for fake_no in xrange(0,len(fits_leftovers)):
    orig_no = fits_name_list.index(fits_leftovers[fake_no])
    bck_sav = os.path.join(backfolder,'soho_' + str(exptimes[orig_no]) + '_background.fits')
    if not os.path.exists(bck_sav):
        if exptimes[orig_no] < 8:
            bck_fits = (exptimes[orig_no]-0.1)/4*(back_4p1 - 390) + 390
        elif exptimes[orig_no] < 19:
            bck_fits = (exptimes[orig_no]-0.1)/4*(back_4p1 - 390) + 390
        else:
            bck_fits = (exptimes[orig_no]-0.1)/4*(back_4p1 - 390) + 390
        fits.writeto(bck_sav, bck_fits, clobber=True)
        
    fits_out_loc = os.path.join(fits_others_save, fits_leftovers[fake_no]) 
    fits_out_bck = os.path.join(fits_others_bgkr, fits_leftovers[fake_no]) 
    fits_in_loc = os.path.join(fits_out_folder, fits_leftovers[fake_no]) 
    data_raw = fits.getdata(fits_in_loc, header=False) 
    fits.writeto(fits_out_loc, data_raw, clobber=True)
    data_bckr = np.clip(data_raw-bck_fits,0,9999999999999999).astype(int)
    fits.writeto(fits_out_bck,data_bckr,clobber=True)
