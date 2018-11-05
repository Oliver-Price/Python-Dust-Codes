# -*- coding: utf-8 -*-
import os
import skimage.transform
import imageio
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
    
nusize = 100

imagedir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour'

smalldir = os.path.join(imagedir,"Small")
if not os.path.exists(smalldir ): os.makedirs(smalldir)

biglist = os.listdir(imagedir)
biglist = [s for s in biglist if "." in s]
bign = len(biglist)

for bid in range(bign):

    image_in = os.path.join(imagedir,biglist[bid])
    im = imageio.imread(image_in)
    
    scalefactor = max(im.shape[0],im.shape[1])/nusize
    
    new_0 = round(im.shape[0]/scalefactor)
    new_1 = round(im.shape[1]/scalefactor)
    
    im_small = (255*skimage.transform.resize(im,(new_0,new_1,3))).astype(int)
    
    if im_small.shape[0] < im_small.shape[1]:
        
       border = im_small.shape[1] - im_small.shape[0]
       
       b = border//2
       t = b + im_small.shape[0]
       
       im_out = np.empty((nusize,nusize,3),dtype=int)
       
       im_out[...,0].fill(int(im_small[...,0].mean()))
       im_out[...,1].fill(int(im_small[...,1].mean()))
       im_out[...,2].fill(int(im_small[...,2].mean()))
       
       im_out[b:t,...] = im_small
    
    if im_small.shape[0] > im_small.shape[1]:
        
       border = im_small.shape[0] - im_small.shape[1]
       
       b = border//2
       t = b + im_small.shape[1]
       
       im_out = np.empty((nusize,nusize,3),dtype=int)
       
       im_out[...,0].fill(int(im_small[...,0].mean()))
       im_out[...,1].fill(int(im_small[...,1].mean()))
       im_out[...,2].fill(int(im_small[...,2].mean()))
       
       im_out[:,b:t,:] = im_small
      
    image_out = os.path.join(smalldir,biglist[bid].split(".")[0]+'.png')
    print(image_out)
    imageio.imwrite(image_out, im_out)
       