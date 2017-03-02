# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

backsave_8p1 = r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\8.1\background_subtracted\20111215_150609_Clear.fits'
backsave_4p1 = r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\4.1\background_subtracted\20111215_170605_Clear.fits'
backsave_19p1 = r'C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\19.1\background_subtracted\20111215_151805_Clear.fits'

values8 = fits.getdata(backsave_8p1, header=False)
values4 = fits.getdata(backsave_4p1, header=False)
values19 = fits.getdata(backsave_19p1, header=False)

x = np.arange(1024)
y8 = values8[512,:]#-values8.min()
y4 = values4[512,:]#-values4.min()
y19 = values19[512,:]#-values19.min()

#y19o4 = y19/y4*4.1/19.1
#y19o8 = y19/y8
#plt.plot(x, y19o8)
#plt.scatter(y4,y4)
#plt.scatter(y4,y8)
#plt.scatter(y4,y19)
#
#y4 = values4[412,:]
#y8 = values8[412,:]
up4 = y4*(19/4)
up8 = y8*(19/8)

#plt.plot(x, y8)
#plt.plot(x, y4)
plt.plot(x, y19)
plt.plot(x, up8)
plt.plot(x, up4)

#plt.plot(x, y19)

plt.show()