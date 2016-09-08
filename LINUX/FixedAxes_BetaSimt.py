#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#with a fixed axis for time over
#*************************************************************

import os
import sys
import time
import pickle
import string
import datetime
import astropy.time
from astropy import wcs
import numpy as np
from astropy.io import fits
import matplotlib.path as mplPath
from PIL import Image, ImageDraw, ImageFont
from dust_functions_LINUX import orb_vector, beta2ypix, simt2xpix, \
radec_slim

#hard coded parameters - can script these in for softer coding
datafolder = '/unsafe/users/op2/panstarrs/stereo_a/'
fpsav = '/unsafe/users/op2/panstarrs/fpsimsav/'
orbitdir = '/unsafe/users/op2/panstarrs/orbitdata/'
pngsav = '/unsafe/users/op2/panstarrs/dustpngs/'
fixaxsav = '/unsafe/users/op2/panstarrs/fixedaxispngs/'
fitssav = '/unsafe/users/op2/panstarrs/dustfits/'
mapsav = '/unsafe/users/op2/panstarrs/pymapsav/'
txtmaps = '/unsafe/users/op2/panstarrs/txtmappings/'
betau = 3.0
betal = 0.02
bno = 600
simtu = 12.0
simtl = 0.5
tno = 600
comdenom = 'C2011L4'
comname = 'Pan-STARRS'
obsloc = 'Stereo B'

#THESE NEED TO BE CHANGED
axlotime = astropy.time.Time(datetime.datetime(2013,3,3,0,0,0))        
axhitime = astropy.time.Time(datetime.datetime(2013,3,11,12,0,0))

#choosing comet data to us
image_list = sorted(os.listdir(mapsav))
image_total = len(image_list)

for image_id in range(170, image_total):

	image_procname = image_list[image_id].split('_mapped')[0]
	image_basename = image_procname[:21]
	image_mapped = os.path.join(mapsav,image_procname + "_mapped.npy")
	image_fixax_img = os.path.join(fixaxsav,image_procname + "_fixedaxisimg.png")	
	image_simressave = os.path.join(fpsav, image_procname + ".npy")
	print (image_basename)
	
	mapping_exists = os.path.exists(image_mapped)
	
	if (mapping_exists == True): #load if it has
		srcolors = np.load(image_mapped)	
		
	elif (mapping_exists == False): 
		print ("No mapped image exists")
		break
	
	simres = np.load(image_simressave)

	#get image time from filename
	csec = int(image_basename[13:15])
	cmin = int(image_basename[11:13])
	chour = int(image_basename[9:11])
	cday = int(image_basename[6:8])
	cmonth = int(image_basename[4:6])
	cyear = int(image_basename[0:4])
	ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
				                               chour , cmin, csec))
		                               
	if not os.path.exists(image_fixax_img):
				                               	    
		t2sfu = ctime.jd - axlotime.jd
		t2sfl = ctime.jd - axhitime.jd
		b1sfu = float('%.1g' % betau)
		b1sfl = float('%.1g' % betal)
    
		pixhi = 600
		pixwt = 1000
		border = 100
		hscle = pixhi/(betau - betal)
		wscle = pixwt/(axhitime.jd - axlotime.jd)

		dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
				                  	pixhi+int(3*border)),(0,0,0,255))
		d = ImageDraw.Draw(dustimg)

		imgmax = srcolors[:,:,0].max()
		if imgmax < 255: imgmax = 255
        
		low = 20000
		hih = 100000
		
		sync_int_list = np.zeros((tno-1), dtype = int)
			
		for tval in range(0,tno-1):
		
			simt_values = np.linspace(simtl, simtu, tno)
			
			sync_int_list[tval] = ((int(np.sum(simres[tval,:,12]) >= 1) + 
												int(np.sum(simres[tval+1,:,12]) >= 1) +
												int(simt_values[tval] <= t2sfu) +
												int(simt_values[tval] >= t2sfl)) == 4)
        
		for ta in np.nonzero(sync_int_list)[0].tolist():
		
			balist1 = np.where(simres[ta+1,:,12] == 1)[0][:-1]
			balist0 = np.where(simres[ta,:,12] == 1)[0][:-1]
			balist = np.intersect1d(balist0,balist1).tolist()
		
			for ba in balist:
			
				fillval = sorted([1, srcolors[ta,ba,0], 9999999999])[1]
				fillco = int(round(255*(np.log10(fillval) - np.log10(low))*
								1 / (np.log10(hih) - np.log10(low))))
				fillco = sorted([0, fillco, 255])[1]
				b1 = beta2ypix(simres[ta,ba,1], border, pixhi, b1sfl, hscle)
				t1 = simt2xpix(simres[ta,ba,0], border, pixwt, t2sfl, wscle)
				b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, b1sfl, hscle)
				t2 = simt2xpix(simres[ta,ba+1,0], border, pixwt, t2sfl, wscle)
				b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, b1sfl, hscle)
				t3 = simt2xpix(simres[ta+1,ba+1,0], border, pixwt, t2sfl, wscle)
				b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, b1sfl, hscle)
				t4 = simt2xpix(simres[ta+1,ba,0], border, pixwt, t2sfl, wscle)
				a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
				,fill=(fillco,fillco,fillco,255))
             
		a = d.polygon([(border,border),(border*2+pixwt,border), \
        (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
        outline = (255,255,255,255))
    
		bmajticks = np.array([0.01,0.1,0.5,1,2,3])
		bminticks = np.array([0.2,0.3,0.4,0.6,0.7,0.8,0.9,1.2,1.4,1.6,1.8,2.2,2.4,2.6,2.8,])
		
		tmajticks = np.arange(t2sfl, t2sfu, 2)
		tminticks = np.arange(t2sfl, t2sfu)
		
		bminticlocs = beta2ypix(bminticks, border, pixhi, b1sfl, hscle)
		bmajticlocs = beta2ypix(bmajticks, border, pixhi, b1sfl, hscle)
		tmajticlocs = simt2xpix(tmajticks, border, pixwt, t2sfl, wscle)
		tminticlocs = simt2xpix(tminticks, border, pixwt, t2sfl, wscle)
 
		majt = 20  #major tick length
		mint = 10  #minor tick length
		xaxis = pixhi + border*2
		fontloc = r'/home/op2/python/lucon.ttf'
		fnt = ImageFont.truetype(fontloc, 20)
		dtime = astropy.time.TimeDelta(1, format='jd')
 
		for div in range(0, (np.size(bminticlocs))): #beta axis minor ticks
			b = d.line([(border+mint,bminticlocs[div]),(border,bminticlocs[div])]
			,fill = (255,255,255,255))
			
		for div in range(0, (np.size(tminticlocs))): #simt axis minor ticks
			b = d.line([(tminticlocs[div],xaxis-mint),(tminticlocs[div],xaxis)]
			,fill = (255,255,255,255))
 
		for div in range(0, (np.size(tmajticlocs))): #simt axis major ticks
			b = d.line([(tmajticlocs[div],xaxis-majt),(tmajticlocs[div],xaxis)],
			fill = (255,255,255,255))
			ticktime = ctime - dtime*tmajticks[div]
			tick = ticktime.isot.replace('T','\n')[0:16]
			d.text((tmajticlocs[div] - len(tick)*3,xaxis + 10), \
			tick, font=fnt, fill=(255,255,255,255))
 
		for div in range(0, (np.size(bmajticlocs))): #beta axis major ticks
			b = d.line([(border+majt,bmajticlocs[div]),(border,bmajticlocs[div])]
			,fill = (255,255,255,255))
			tick = str(bmajticks[div])
			d.text((border - len(tick)*5 - 40,bmajticlocs[div] - 10 ), \
			tick, font=fnt, fill=(255,255,255,255))
     
		#axis labels
		d.text((1.5*border + pixwt*0.5 - 150,pixhi + 2.7*border), \
		"Date/Time of Ejection", font=fnt, fill=(255,255,255,255))
		d.text((0.25*border,0.75*border), \
		"Beta", font=fnt, fill=(255,255,255,255))

		#plot title
		plttitle = (comdenom + ' ' + comname + ' from ' + 
		obsloc + ' from date: ' + ctime.isot[0:16].replace('T',' at '))
		tfnt = ImageFont.truetype(fontloc, 30)
		d.text((1.5*border + pixwt*0.5 - len(plttitle)*8 - 90,.35*border), \
		plttitle, font=tfnt, fill=(255,255,255,255))
		
		if (t2sfl <= 0 <= t2sfu):
			image_line_x = simt2xpix(0, border, pixwt, t2sfl, wscle)
			b = d.line([(image_line_x,border),(image_line_x,2*border+pixhi)]
			,fill = (255,0,0,255))
			d.text((image_line_x - 160,1.25*border), \
			"Time of Image", font=fnt, fill=(255,0,0,255))


		#dustimg.show()
		dustimg.save(image_fixax_img,'png')

