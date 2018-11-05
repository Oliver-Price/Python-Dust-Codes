# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import os
import pandas as pd

xmldir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\Labels'
xlist = os.listdir(xmldir)

d = {'file':[],'Nuc X':[],'Nuc Y':[],'Ion End X':[],'Ion End Y':[]}
df = pd.DataFrame(data = d)

for f in range(len(xlist)):
    
    exml = os.path.join(xmldir,xlist[f])
    
    tree = ET.parse(exml)
    root = tree.getroot()
    save = ''
    
    for child in root:
        for grandchild in child:
            if grandchild.tag == 'name':
                if grandchild.text == 'nucleus':
                    save = 'n'
                elif grandchild.text == 'ion tail':
                    save = 'i'
                else:
                    save = ''
            if grandchild.tag == 'bndbox':
                if save == 'n':
                    yNmin = float(grandchild[0].text)
                    xNmin = float(grandchild[1].text)
                    yNmax = float(grandchild[2].text)
                    xNmax = float(grandchild[3].text)
                if save == 'i':
                    yImin = float(grandchild[0].text)
                    xImin = float(grandchild[1].text)
                    yImax = float(grandchild[2].text)
                    xImax = float(grandchild[3].text)
                    
    yNave = 0.5*(yNmin + yNmax)
    xNave = 0.5*(xNmin + xNmax)
    
    if abs(yImin - yNave) < abs(yImax - yNave):
        Nuc_y = yImin
        Ion_y = yImax
    else:
        Nuc_y = yImax
        Ion_y = yImin
        
    if abs(xImin - xNave) < abs(xImax - xNave):
        Nuc_x = xImin
        Ion_x = xImax
    else:
        Nuc_x = xImax
        Ion_x = xImin

    fid = xlist[f].split('.')[0]
        
    df_line = pd.DataFrame(data = {'file':[fid],'Nuc X':[Nuc_x],'Nuc Y':[Nuc_y],'Ion End X':[Ion_x],'Ion End Y':[Ion_y]})
    df = pd.concat((df,df_line))

df.to_pickle('C:\PhD\Python\Save_Data\df.pickle')