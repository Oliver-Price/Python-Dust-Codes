# -*- coding: utf-8 -*-
#Original:
#gamma=4.5,h=0.4,k=1.4
#gvals = [3.5,4.5,5.5]
#hvals = [0.2,0.3,0.4]
#kvals = [1.4,1.6,1.8]

import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
from scikit_enhance import mgn
import imageio

def rgb2gray(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

fitsdir = r'C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Earth'
dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

outdir = r'C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Earth\mgntest'

gvals = [3.2,4.2,5.2]
hvals = [0.2,0.4,0.6]
kvals = [1.0,1.5,2.0]

fits_id = 2
fitstemp = os.path.join(fitsdir, fits_list[fits_id])
hdulist = fits.open(fitstemp)
data = (hdulist[0].data)
img_size = np.shape(data)
data = np.swapaxes(data,0,2)
data = rgb2gray(data)
data = np.nan_to_num(data).astype(float)
for g in gvals:
    for h in hvals:
        for k in kvals:
            outname = "gamma_" + str(g) + "_h_" + str(h) + "_k_" + str(k) + ".png"
            enhanced = mgn(data,k=k,gamma=g,h=h)
            outfile = os.path.join(outdir,outname)
            imageio.imwrite(outfile, enhanced)




