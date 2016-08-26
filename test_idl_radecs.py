#testing idl soho radec
from astropy.io import fits
import numpy as np
from plot_functions import ra2xpix, dec2ypix, setaxisup
from PIL import Image, ImageDraw, ImageFont

fitsin = r"C:\Users\op2\Downloads\01671892\35170060.fts"

hdulist = fits.open(fitsin)
fits_vals = (hdulist[0].data)

raloc = r"C:\PhD\IDL\IDL_codes\Testcodes\ra.txt"
decloc = r"C:\PhD\IDL\IDL_codes\Testcodes\dec.txt"

ra = np.genfromtxt(raloc, delimiter = ",")
dec = np.genfromtxt(decloc, delimiter = ",")

ramin = np.amin(ra)
ramax = np.amax(ra)
decmin = np.amin(dec)
decmax = np.amax(dec)

if ramax < 0:
    ra_m = np.copy(ra) + 360
    rafmin = np.amin(ra_m)
    rafmax = np.amax(ra_m)
elif (ramax-ramin) > 270:
    ra_m = np.copy(ra)
    circlocs = np.where(ra_m < 180)
    ra_m[circlocs] = ra_m[circlocs] + 360
    rafmin = np.amin(ra_m)
    rafmax = np.amax(ra_m)
else:
    ra_m = ra
    rafmin = ramin
    rafmax = ramax

backgr_fill = (0,0,0,255)
featur_fill = (255,255,255,255)
pixheight = 800
pixwidth = int(pixheight*(rafmax - rafmin)/(decmax - decmin))
border = 100
scale = pixheight/(decmax - decmin)
imgwidth = pixwidth+int(4*border)
imgheight = pixheight+int(3*border)
comimg = Image.new('RGBA', ( imgwidth , imgheight ) ,backgr_fill)
d = ImageDraw.Draw(comimg)

xa = fits_vals.shape[0]
ya = fits_vals.shape[1]
low = 2e-10
hih = 2e-09

for x in xrange(0, ya-2):
    for y in xrange(0, xa-2):
        fillco = int(round(255*(np.log10(fits_vals[x,y]) - np.log10(low))*
								1 / (np.log10(hih) - np.log10(low))))
        fillco = sorted([0, fillco, 255])[1]
        ra1 = ra2xpix(ra_m[x,y],border,pixwidth,rafmin,scale)
        ra2 = ra2xpix(ra_m[x+1,y],border,pixwidth,rafmin,scale)
        ra3 = ra2xpix(ra_m[x+1,y+1],border,pixwidth,rafmin,scale)
        ra4 = ra2xpix(ra_m[x,y+1],border,pixwidth,rafmin,scale)
        dec1 = dec2ypix(dec[x,y],border,pixheight,decmin,scale)
        dec2 = dec2ypix(dec[x+1,y],border,pixheight,decmin,scale)
        dec3 = dec2ypix(dec[x+1,y+1],border,pixheight,decmin,scale)
        dec4 = dec2ypix(dec[x,y+1],border,pixheight,decmin,scale)
        a = d.polygon([(ra1,dec1),(ra2,dec2),(ra3,dec3),(ra4,dec4)] ,\
        fill=(fillco,fillco,fillco,255))

a = d.polygon([(border,border),(border*2+pixwidth,border), \
(border*2+pixwidth,border*2+pixheight),(border,border*2+pixheight)], \
outline = featur_fill)

axisdata = setaxisup(rafmax,rafmin,decmax,decmin,border,pixheight,pixwidth,scale)

majt = 20  #major tick length
mint = 10  #minor tick length
rdaxis = pixheight + border*2

fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
fnt = ImageFont.truetype(fontloc, 20)

for div in xrange(0, (np.size(axisdata[1]))): #RA axis major ticks
    b = d.line([(axisdata[1][div],rdaxis-majt),(axisdata[1][div],rdaxis)],\
    fill = featur_fill)
    tick = str(axisdata[0][div]%360)
    d.text((axisdata[1][div] - len(tick)*5,rdaxis + 10), \
    tick, font=fnt, fill= featur_fill)
    
for div in xrange(0, (np.size(axisdata[2]))): #RA axis minor ticks
    b = d.line([(axisdata[2][div],rdaxis-mint),(axisdata[2][div],rdaxis)],\
    fill= featur_fill)

for div in xrange(0, (np.size(axisdata[4]))): #DEC axis major ticks
    b = d.line([(border+majt,axisdata[4][div]),(border,axisdata[4][div])],\
    fill= featur_fill)
    tick = str(axisdata[3][div])
    d.text((border - len(tick)*5 - 40,axisdata[4][div] - 10 ), \
    tick, font=fnt, fill=(255,255,255,128))
    
for div in xrange(0, (np.size(axisdata[5]))): #DEC axis minor ticks
    b = d.line([(border+mint,axisdata[5][div]),(border,axisdata[5][div])],\
    fill= featur_fill)

#axis labels
d.text((1.5*border + pixwidth*0.5 - 145,pixheight + 2.5*border), \
"Right Ascension (Degrees)", font=fnt, fill= featur_fill)
d.text((0.25*border - 10,0.75*border - 20), \
"Declination (Degrees)", font=fnt, fill= featur_fill)

comimg.show()