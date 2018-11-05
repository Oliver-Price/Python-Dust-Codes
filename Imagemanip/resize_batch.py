# -*- coding: utf-8 -*-

from PIL import Image
from resizeimage import resizeimage
import os

file_dir = r'C:\PhD\Comet_data\Comet_Hale_Bopp\all_images'

im_out_dir = os.path.join(file_dir,'small')
if not os.path.exists(im_out_dir):
    os.makedirs(im_out_dir)

imgs = os.listdir(file_dir)

for i in range(len(imgs)):
    imagepath = os.path.join(file_dir,imgs[i])
    
    with open(imagepath, 'r+b') as img:
        with Image.open(img) as image:
            re_img  = resizeimage.resize_thumbnail(image, [100, 100])
    
    '''       
    try:
        im_out_base = imgs[i].split('_')[-1]
    except:
        im_out_base = imgs[i]
    '''     
    
    im_out_base = 'HaleBopp' + str(i) + '.png'

    im_out = os.path.join(im_out_dir,im_out_base)
        
    re_img.save(im_out, re_img.format)