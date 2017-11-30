import os
from astropy.io import fits
import numpy as np
from PIL import Image, ImageDraw, ImageFont

fitsdir = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1-diff'
pngdir = os.path.join(fitsdir,'pngs2')
if not os.path.exists(pngdir): os.makedirs(pngdir)
#mloc = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\processing\filter_MASK.png'
#maskimg = Image.open(mloc)
#maska = np.flipud(np.array(maskimg)[:,:,0]/255)

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

low = -150; hih=150
methodlog = False

fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
fnt = ImageFont.truetype(fontloc, 50)
featur_fill = (255,0,0,255)

for fits_id in range(40,fits_total):
    
        fitstemp = os.path.join(fitsdir, fits_list[fits_id])
        pngtemp = os.path.join(pngdir, (fits_list[fits_id].split(".")[0] +
        "_" + str(low) + "_" + str(hih) + ".png"))
        
        if not os.path.exists(pngtemp):
            
            #OP filename method
            cmin = os.path.basename(fitstemp)[11:13]
            chour = os.path.basename(fitstemp)[9:11]
            cday = os.path.basename(fitstemp)[6:8]
            cmonth = os.path.basename(fitstemp)[4:6]
            cyear = os.path.basename(fitstemp)[0:4]
  
            #YR filename method          
#            cmin =os.path.basename(fitstemp)[21:23]
#            chour = os.path.basename(fitstemp)[19:21]
#            cday = os.path.basename(fitstemp)[16:18]
#            cmonth = os.path.basename(fitstemp)[13:15]
#            cyear = os.path.basename(fitstemp)[8:12]

            hdulist = fits.open(fitstemp)
            data = (hdulist[0].data)
            img_size = np.shape(data)
            data = np.nan_to_num(data)#*maska
            
            pngimg = Image.new('RGBA', (img_size[0], img_size[1] ) , (0,0,0,255) )
            d = ImageDraw.Draw(pngimg)
            
            
            if methodlog == True:
                grad = np.log10(hih/low)                        
                filldat = (np.clip(255.0/grad*np.log10(data/low),0,255)).astype(int)
                for x in range(0,img_size[0]):
                    for y in range(0,img_size[1]):          
                        d.point((x,y), fill = (filldat[x,y],filldat[x,y],filldat[x,y]) )
                  
            if methodlog == False:
                grad = (hih - low)                     
                filldat = (np.clip(255.0/grad*(data-low),0,255)).astype(int)
                for x in range(0,img_size[0]):
                    for y in range(0,img_size[1]):          
                        d.point((x,y), fill = (filldat[x,y],filldat[x,y],filldat[x,y]) )

            times = (cday + '/' + cmonth + '/' + cyear + '\n' + chour + ':' + cmin)
            
            #d.text((650,650), times , font=fnt, fill= featur_fill)
            pngimg.save(pngtemp,'png')
        print (fits_id)