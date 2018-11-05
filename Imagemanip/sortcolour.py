# -*- coding: utf-8 -*-
import os
import numpy as np
import scipy.misc

imagedir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images'

bwdir = os.path.join(imagedir,"Black and White")
if not os.path.exists(bwdir): os.makedirs(bwdir)
rgbdir = os.path.join(imagedir,"Colour")
if not os.path.exists(rgbdir): os.makedirs(rgbdir)

biglist = os.listdir(imagedir)
biglist = [s for s in biglist if "." in s]
bign = len(biglist)

for x in range(bign):
    image_in = os.path.join(imagedir,biglist[x])
    image = scipy.misc.imread(image_in)
    if image.ndim == 3:
        if np.all(image[...,0] == image[...,1]) == False:
            image_out = os.path.join(rgbdir,biglist[x])            
        else:
            image_out = os.path.join(bwdir,biglist[x])
    elif image.ndim == 2:
        image_out = os.path.join(bwdir,biglist[x])
    else:
        pass
    
    os.rename(image_in,image_out)