# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import os

xmldir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\Labels'
xlist = os.listdir(xmldir)

for f in range(len(xlist)):
    
    nuclabelled = False
    ionlabelled = False
    
    exml = os.path.join(xmldir,xlist[f])
    
    tree = ET.parse(exml)
    root = tree.getroot()
        
    for child in root:
        for grandchild in child:
            if grandchild.tag == 'name':
                nuclabelled = (grandchild.text == 'nucleus') or nuclabelled
                ionlabelled = (grandchild.text == 'ion tail') or ionlabelled
    
    if not (nuclabelled and ionlabelled):
        print (exml)