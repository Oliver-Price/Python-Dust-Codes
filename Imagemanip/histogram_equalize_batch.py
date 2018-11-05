# -*- coding: utf-8 -*-

from PIL import Image
from resizeimage import resizeimage
import os
import sys
import numpy as np
from skimage import exposure
from skimage.color.adapt_rgb import adapt_rgb, each_channel, hsv_value
import scipy.misc

@adapt_rgb(hsv_value)
def equal_hsv(image):
    return exposure.equalize_hist(image)

file_dir = r'C:\PhD\Comet_data\Comet_Hale_Bopp\all_images'

im_out_dir = os.path.join(file_dir,'equalised')
if not os.path.exists(im_out_dir):
    os.makedirs(im_out_dir)

imgs = os.listdir(file_dir)
imgs = [x for x in imgs if "." in x]

for i in range(len(imgs)):
    imagepath = os.path.join(file_dir,imgs[i])
    
    with open(imagepath, 'r+b') as img:
        with Image.open(img) as image:
            re_img  = resizeimage.resize_thumbnail(image, [100, 100])
    
    pix = np.array(re_img)
    
    img_eq = equal_hsv(pix)

    im_out_base = 'HaleBopp' + str(i) + '.png'

    im_out = os.path.join(im_out_dir,im_out_base)
     
    scipy.misc.imsave(im_out,img_eq)