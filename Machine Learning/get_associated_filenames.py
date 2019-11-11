# -*- coding: utf-8 -*-
import os

fitsdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\fits'
origdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

#for fits_id in range(0,fits_total):