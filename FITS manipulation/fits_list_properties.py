from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt

fits_wcs_folder = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue\processing\OLD\small1"
fits_wcs_list = os.listdir(fits_wcs_folder)
fits_wcs_list = [s for s in fits_wcs_list if ".f" in s]

exptimes = np.zeros((len(fits_wcs_list)),dtype=float)
maxval = np.zeros((len(fits_wcs_list)),dtype=float)
obstime = []

for fits_no in range(0,len(fits_wcs_list)):
    
    fits_wcs_loc = os.path.join(fits_wcs_folder, fits_wcs_list[fits_no])
    
    data_wcs, header_wcs = fits.getdata(fits_wcs_loc, header=True)

    exptimes[fits_no] = header_wcs['EXPTIME']
    maxval[fits_no] = np.max(data_wcs)
    obstime .append(header_wcs['DATE-OBS'])     
    
#plt.scatter(exptimes,maxval)
#plt.yscale('log')
#plt.xscale('log')

