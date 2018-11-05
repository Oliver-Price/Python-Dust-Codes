#PLOT MULTIPLE IMAGE WITH SAME RA AND DEC AXIS

import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
from PIL import Image, ImageDraw, ImageFont

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from FP_plot_functions import setaxisup, plotpixel2, ra2xpix, dec2ypix
from conversion_routines import fixwraps

fits1 = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1\20070115_120100_s4h1A.fts'
fits2 = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-2\20070115_120100_s4h2A.fts'
fits3 = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\20070115_115313_Clear.fits'

onedimg1 = fits.open(fits1)
try:
    w1 = wcs.WCS(onedimg1[0].header, key = 'A', naxis=2)
except:
    w1 = wcs.WCS(onedimg1[0].header, naxis=2)
colours1 = np.nan_to_num(onedimg1[0].data)
        
#make a 2xN array of all pixel locations
ya1 = colours1.shape[0]
xa1 = colours1.shape[1]
xv1, yv1= np.meshgrid(np.arange(xa1), np.arange(ya1))
xv1 = np.reshape(xv1, (xa1*ya1,1))
yv1 = np.reshape(yv1, (xa1*ya1,1))
coords1 = np.zeros((xa1*ya1,2))
coords1[:,0] = xv1[:,0]
coords1[:,1] = yv1[:,0]

#convert each pixel locaion to an RA and DEC array
radecs1 = w1.wcs_pix2world(coords1, 0)
ra1 = np.reshape(radecs1[:,0], (ya1,xa1))
dec1 = np.reshape(radecs1[:,1], (ya1,xa1))

#find minimum/maximum values
ramin1 = np.amin(ra1)
ramax1 = np.amax(ra1)
decmin1 = np.amin(dec1)
decmax1 = np.amax(dec1)
        
[ra_m1, rafmin1, rafmax1, bool_val1] = fixwraps(ra1, ramax1, ramin1)

onedimg2 = fits.open(fits2)
try:
    w2 = wcs.WCS(onedimg2[0].header, key = 'A', naxis=2)
except:
    w2 = wcs.WCS(onedimg2[0].header, naxis=2)
colours2 = np.nan_to_num(onedimg2[0].data)
        
#make a 2xN array of all pixel locations
ya2 = colours2.shape[0]
xa2 = colours2.shape[1]
xv2, yv2= np.meshgrid(np.arange(xa2), np.arange(ya2))
xv2 = np.reshape(xv2, (xa2*ya2,1))
yv2 = np.reshape(yv2, (xa2*ya2,1))
coords2 = np.zeros((xa2*ya2,2))
coords2[:,0] = xv2[:,0]
coords2[:,1] = yv2[:,0]

#convert each pixel locaion to an RA and DEC array
radecs2 = w2.wcs_pix2world(coords2, 0)
ra2 = np.reshape(radecs2[:,0], (ya2,xa2))
dec2 = np.reshape(radecs2[:,1], (ya2,xa2))

#find minimum/maximum values
ramin2 = np.amin(ra2)
ramax2 = np.amax(ra2)
decmin2 = np.amin(dec2)
decmax2 = np.amax(dec2)
        
[ra_m2, rafmin2, rafmax2, bool_val2] = fixwraps(ra2, ramax2, ramin2)

onedimg3 = fits.open(fits3)
try:
    w3 = wcs.WCS(onedimg3[0].header, key = 'A', naxis=2)
except:
    w3 = wcs.WCS(onedimg3[0].header, naxis=2)
colours3 = np.nan_to_num(onedimg3[0].data)
        
#make a 3xN array of all pixel locations
ya3 = colours3.shape[0]
xa3 = colours3.shape[1]
xv3, yv3= np.meshgrid(np.arange(xa3), np.arange(ya3))
xv3 = np.reshape(xv3, (xa3*ya3,1))
yv3 = np.reshape(yv3, (xa3*ya3,1))
coords3 = np.zeros((xa3*ya3,2))
coords3[:,0] = xv3[:,0]
coords3[:,1] = yv3[:,0]

#convert each pixel locaion to an RA and DEC array
radecs3 = w3.wcs_pix2world(coords3, 0)
ra3 = np.reshape(radecs3[:,0], (ya3,xa3))
dec3 = np.reshape(radecs3[:,1], (ya3,xa3))

#find minimum/maximum values
ramin3 = np.amin(ra3)
ramax3 = np.amax(ra3)
decmin3 = np.amin(dec3)
decmax3 = np.amax(dec3)
        
[ra_m3, rafmin3, rafmax3, bool_val3] = fixwraps(ra3, ramax3, ramin3)

featur_fill = (0,0,0,255)    
backgr_fill = (255,255,255,255)

trafmin = sorted([rafmin1,rafmin2,rafmin3])[0]
trafmax = sorted([rafmax1,rafmax2,rafmax3])[2]
tdecmin = sorted([decmin1,decmin2,decmin3])[0]
tdecmax = sorted([decmax1,decmax2,decmax3])[2]

pixheight = 800
pixwidth = int(pixheight*(trafmax - trafmin)/(tdecmax - tdecmin))
border = 100
scale = pixheight/(tdecmax - tdecmin)
imgwidth = pixwidth+int(4*border)
imgheight = pixheight+int(3*border)
comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,backgr_fill)

d = ImageDraw.Draw(comimg)
low1 = 10000; hih1 = 1500000
low2 = 10000; hih2 = 1500000
low3 = 4.8e-10; hih3 = 1.4e-9        
         
grad1 = (255 / (np.log10(hih1) - np.log10(low1)))   
grad2 = (255 / (np.log10(hih2) - np.log10(low2)))
grad3 = (255 / (np.log10(hih3) - np.log10(low3)))   

colplot1 = (np.clip(np.round((np.log10(np.clip(colours1,1e-20,999999999))- np.log10(low1))*grad1),0,255)).astype(int)
colplot2 = (np.clip(np.round((np.log10(np.clip(colours2,1e-20,999999999))- np.log10(low2))*grad2),0,255)).astype(int)
colplot3 = (np.clip(np.round((np.log10(np.clip(colours3,1e-20,999999999))- np.log10(low3))*grad3),0,255)).astype(int)

imagemask1 = np.ones_like(colours1)[:-1,:-1]
#imagemask2 = np.ones_like(colours2)[:-1,:-1]
imagemask3 = np.ones_like(colours3)[:-1,:-1]
imagemask3[np.where(colours3 == 0)] = 0

mask2loc = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\hi2A_mask.fts'
maskfits = fits.open(mask2loc)
imagemask2 = (maskfits[0].data[::2,::2])[:-1,:-1]

non_zeros_01 = np.where(imagemask1 == 1)[0]
non_zeros_11 = np.where(imagemask1 == 1)[1]
no_points1 = non_zeros_01.size 

non_zeros_02 = np.where(imagemask2 == 1)[0]
non_zeros_12 = np.where(imagemask2 == 1)[1]
no_points2 = non_zeros_02.size 

non_zeros_03 = np.where(imagemask3 == 1)[0]
non_zeros_13 = np.where(imagemask3 == 1)[1]
no_points3 = non_zeros_03.size 

rp1 = ra2xpix(ra_m1, border, pixwidth, trafmin, scale)
dp1 = dec2ypix(dec1, border, pixheight, tdecmin, scale)

rp2 = ra2xpix(ra_m2, border, pixwidth, trafmin, scale)
dp2 = dec2ypix(dec2, border, pixheight, tdecmin, scale)

rp3 = ra2xpix(ra_m3, border, pixwidth, trafmin, scale)
dp3 = dec2ypix(dec3, border, pixheight, tdecmin, scale)    

for xp in range(0,no_points3):
    x = non_zeros_03[xp]
    y = non_zeros_13[xp]
    plotpixel2(d,x,y,rp3,dp3,colplot3[x,y],colplot3[x,y],colplot3[x,y])

for xp in range(0,no_points1):
    x = non_zeros_01[xp]
    y = non_zeros_11[xp]
    plotpixel2(d,x,y,rp1,dp1,colplot1[x,y],colplot1[x,y],colplot1[x,y])
        
for xp in range(0,no_points2):
    x = non_zeros_02[xp]
    y = non_zeros_12[xp]
    plotpixel2(d,x,y,rp2,dp2,colplot2[x,y],colplot2[x,y],colplot2[x,y])

#draws a border       
a = d.polygon([(border,border),(border*2+pixwidth,border), \
(border*2+pixwidth,border*2+pixheight),(border,border*2+pixheight)], \
outline = featur_fill)
        
#most of the dirty stuff is bunged into this function
axisdata = setaxisup(trafmax,trafmin,tdecmax,tdecmin,border,pixheight,pixwidth,scale)

majt = 20  #major tick length
mint = 10  #minor tick length
rdaxis = pixheight + border*2

fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
fnt = ImageFont.truetype(fontloc, 20)
smallfnt = ImageFont.truetype(fontloc, 10)
largefnt = ImageFont.truetype(fontloc, 30)

for div in range(0, (np.size(axisdata[1]))): #RA axis major ticks
    b = d.line([(axisdata[1][div],rdaxis-majt),(axisdata[1][div],rdaxis)],\
    fill = featur_fill)
    tick = str(axisdata[0][div]%360)
    d.text((axisdata[1][div] - len(tick)*5,rdaxis + 10), \
    tick, font=fnt, fill= featur_fill)
    
for div in range(0, (np.size(axisdata[2]))): #RA axis minor ticks
    b = d.line([(axisdata[2][div],rdaxis-mint),(axisdata[2][div],rdaxis)],\
    fill= featur_fill)

for div in range(0, (np.size(axisdata[4]))): #DEC axis major ticks
    b = d.line([(border+majt,axisdata[4][div]),(border,axisdata[4][div])],\
    fill= featur_fill)
    tick = str(axisdata[3][div])
    d.text((border - len(tick)*5 - 40,axisdata[4][div] - 10 ), \
    tick, font=fnt, fill=featur_fill)
    
for div in range(0, (np.size(axisdata[5]))): #DEC axis minor ticks
    b = d.line([(border+mint,axisdata[5][div]),(border,axisdata[5][div])],\
    fill= featur_fill)
    
#axis labels
d.text((1.5*border + pixwidth*0.5 - 145,pixheight + 2.5*border), \
"Right Ascension (Degrees)", font=fnt, fill= featur_fill)
d.text((0.25*border - 10,0.75*border - 20), \
"Declination (Degrees)", font=fnt, fill= featur_fill)

comimg.show()
comimg.save(r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\15th_midday_obs.png','png')