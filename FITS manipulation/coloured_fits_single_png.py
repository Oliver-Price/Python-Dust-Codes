import os
from astropy.io import fits
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib

fitsfile = r"C:\Users\Ollie\Downloads\MLSO\20030218_020207.fts"

hih = 600
low = -1000

hdulist = fits.open(fitsfile)
data = np.rot90((hdulist[0].data),-1)
img_size = np.shape(data)

grad = (hih - low)                     
filldat = (np.clip(1.0/grad*(data-low),0,1))

cmap = matplotlib.cm.get_cmap('hot')
coldat = (cmap(filldat)*255).astype(int)

pngimg = Image.new('RGBA', (img_size[0], img_size[1] ) , (0,0,0,255) )
d = ImageDraw.Draw(pngimg)

for (x, y), element in np.ndenumerate(data):    
    d.point((x,y), fill = (coldat[x,y,0],coldat[x,y,1],coldat[x,y,2]) )
    
pngimg.save(r"C:\PhD\Python\file_directory\MLSO.png",'png')

    