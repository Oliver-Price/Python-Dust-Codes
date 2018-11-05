# -*- coding: utf-8 -*-
from PIL import Image, ImageDraw
import os
import pandas as pd

df = pd.read_pickle('C:\PhD\Python\Save_Data\df.pickle')

file_dir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour'
imgs = os.listdir(file_dir)
imgs = [a for a in imgs if '.' in a]

im_out_dir = os.path.join(file_dir,'labelled')
if not os.path.exists(im_out_dir):
    os.makedirs(im_out_dir)

im_nolbl = os.path.join(file_dir,'nolabel')
if not os.path.exists(im_nolbl):
    os.makedirs(im_nolbl)

for i in range(len(imgs)):
    
    image_out = os.path.join(im_out_dir,imgs[i])
    
    if not os.path.exists(image_out):
    
        imagepath = os.path.join(file_dir,imgs[i])
        
        filebase = imgs[i].split('.')[0]
        
        dfa = df[df['file'].str.match(filebase)]
        
        if dfa.empty == True:
            
            nolblpath = os.path.join(im_nolbl,imgs[i])
            
            os.rename(imagepath,nolblpath)
        
        else:
            
            dfs = dfa.iloc[0]
            
            im = Image.open(imagepath)
            
            d = ImageDraw.Draw(im)       
            
            d.line((dfs[2], dfs[1], dfs[4], dfs[3]), fill=(255,0,0,255))
            
            im.save(image_out)