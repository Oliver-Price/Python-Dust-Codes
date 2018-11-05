# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import os
from PIL import Image, ImageDraw

xmldir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\Labels'
xlist = os.listdir(xmldir)

file_dir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour'
imgs = os.listdir(file_dir)
imgs = [a for a in imgs if '.' in a]

im_out_dir = os.path.join(file_dir,'notext')
if not os.path.exists(im_out_dir):
    os.makedirs(im_out_dir)

for f in range(0,len(xlist)):
    
    exml = os.path.join(xmldir,xlist[f])
    filebase = xlist[f].split('.')[0]
    imgbase = [a for a in imgs if filebase in a][0]
    image_out = os.path.join(im_out_dir,imgbase)
    
    if not os.path.exists(image_out):
        tree = ET.parse(exml)
        root = tree.getroot()
        save = False
        textlist = []    
        
        for child in root:
            for grandchild in child:
                if grandchild.tag == 'name':
                    if grandchild.text == 'text':
                        save = True    
                if grandchild.tag == 'bndbox':
                    if save == True:
                        ymin = float(grandchild[0].text)
                        xmin = float(grandchild[1].text)
                        ymax = float(grandchild[2].text)
                        xmax = float(grandchild[3].text)
                        textlist.append((ymin,ymax,xmin,xmax))
                        save = False
        
        imagepath = os.path.join(file_dir,imgbase)    
        im = Image.open(imagepath)
        d = ImageDraw.Draw(im)
        
        for i in range(len(textlist)):
            (ymin,ymax,xmin,xmax) = textlist[i]
            d.polygon([(ymin,xmin),(ymin,xmax),(ymax,xmax),(ymax,xmin)] ,fill=(0,0,0,255))
            
        im.save(image_out)
        