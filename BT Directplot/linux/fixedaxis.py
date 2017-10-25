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
sys.path.append(r"/home/op2/python/mcnaught/")
from dust_functions_LINUX import orb_obs, logbeta2ypix, simt2xpix, round_to_base, RoundToSigFigs

#hard coded parameters - can script these in for softer coding
fpsav = '/unsafe/users/op2/mcnaught_stereo/HI-1/fpsimsav/narrow3/'
fixaxsav = '/unsafe/users/op2/mcnaught_stereo/HI-1-MGN/fixedaxispngs/narrow3/'
mapsav = '/unsafe/users/op2/mcnaught_stereo/HI-1-MGN/pymapsav/narrow3/'
betau = 3.5
betal = 0.2
bno = 1000
simtu = 16.0
simtl = 2.0
tno = 1000
comdenom = 'C2006P1'
comname = 'McNaught'
obsloc = 'Stereo A' 

axlotime = astropy.time.Time(datetime.datetime(2006,12,27,0,1,0))
axhitime = astropy.time.Time(datetime.datetime(2007,1,14,0,1,0))

#choosing comet data to us
image_list = sorted(os.listdir(mapsav))
image_total = len(image_list)

for image_id in range(%START%,%FINISH%):

	image_procname = image_list[image_id].split('.')[0][:-7]
	image_basename = image_procname[:21] + '_diff'
	image_mapped = os.path.join(mapsav,image_procname + "_mapped.npy")
	image_simressave = os.path.join(fpsav, image_procname+ ".npy")
	image_fixax_img = os.path.join(fixaxsav,image_procname + "_fixedaxisimg.png")	
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
				                               	    
		tuppr = ctime.jd - axlotime.jd
		tlowr = ctime.jd - axhitime.jd
		bl1sf = RoundToSigFigs(betal,1); bu1sf = RoundToSigFigs(betau,1)
	
		pixhi = 600
		pixwt = 1000
		border = 100
		hscle = pixhi/np.log(betau/betal)
		wscle = pixwt/(axhitime.jd - axlotime.jd)

		dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
				                  	pixhi+int(3*border)),(0,0,0,255))
		d = ImageDraw.Draw(dustimg)

		srcolors_slim = np.copy(srcolors[:,:,0])
		srcolors_slim[np.where(simres[:,:,0] >= tuppr)] = 0
		srcolors_slim[np.where(simres[:,:,0] <= tlowr)] = 0

		good_bool = srcolors_slim[1:,1:] + srcolors_slim[:-1,:-1] + srcolors_slim[1:,:-1] + srcolors_slim[:-1,1:]
		good_locs = np.where(good_bool==4)

		low = -0.7
		hih = 1.25
		sr_fill = np.clip(255.0/(hih-low)*(srcolors[:,:,1]-low),0,255).astype(int)
        
		for x in range(0, np.size(good_locs[0])):
			ta = good_locs[0][x]; ba = good_locs[1][x]
			b1 = logbeta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
			t1 = simt2xpix(simres[ta,ba,0], border, pixwt, tlowr, wscle)
			b2 = logbeta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
			t2 = simt2xpix(simres[ta,ba+1,0], border, pixwt, tlowr, wscle)
			b3 = logbeta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
			t3 = simt2xpix(simres[ta+1,ba+1,0], border, pixwt, tlowr, wscle)
			b4 = logbeta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
			t4 = simt2xpix(simres[ta+1,ba,0], border, pixwt, tlowr, wscle)
			a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
			,fill=(sr_fill[ta,ba],sr_fill[ta,ba],sr_fill[ta,ba],255))
             
		a = d.polygon([(border,border),(border*2+pixwt,border), \
        (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
        outline = (255,255,255,255))
    
		tdivmajors = np.array([1.,2.,3.,4.,5.])
		tdivnos = (1/tdivmajors) * (tuppr - tlowr)
		tnodivs = 6
		tdividx = np.where((tdivnos <= tnodivs)==True)[0][0]
		tmajticks = np.arange(tuppr, tlowr-0.0001, -tdivmajors[tdividx])

		tdivminors = np.array([0.5,1,1,1,1])
		tminticks = np.arange(tuppr, tlowr-0.0001, -tdivminors[tdividx])
		tminticks = np.setdiff1d(tminticks,tmajticks)
		
		bdivmajors = np.array([0.1,0.2,0.5,1,2])
		bdivnos = (1/bdivmajors) * (bu1sf - bl1sf)
		bnodivs = 10
		bdividx = np.where((bdivnos <= bnodivs)==True)[0][0]
		bu2majdv = round_to_base(betau, bdivmajors[bdividx])
		bl2majdv = round_to_base(betal, bdivmajors[bdividx])   
		bmajticks = np.arange(bu2majdv, bl2majdv-0.0001, -bdivmajors[bdividx])
		bmajticks = bmajticks[bmajticks >= betal-0.0002]
		bmajticks = np.round(bmajticks,2)

		bdivminors = np.array([0.02,0.05,0.1,0.2,0.5])
		bu2mindv = round_to_base(betau, bdivminors[bdividx])
		bl2mindv = round_to_base(betal, bdivminors[bdividx])   
		bminticks = np.arange(bu2mindv, bl2mindv-0.0001, -bdivminors[bdividx])
		bminticks = np.setdiff1d(bminticks,bmajticks)
		bminticks = bminticks[bminticks >= betal-0.0002]
		bminticks = np.round(bminticks,2)

		bminticlocs = logbeta2ypix(bminticks, border, pixhi, betal, hscle)
		bmajticlocs = logbeta2ypix(bmajticks, border, pixhi, betal, hscle)
		tmajticlocs = simt2xpix(tmajticks, border, pixwt, tlowr, wscle)
		tminticlocs = simt2xpix(tminticks, border, pixwt, tlowr, wscle)
 
		majt = 20  #major tick length
		mint = 10  #minor tick length
		xaxis = pixhi + border*2
		fontloc = r'/home/op2/python/lucon.ttf'
		fnt = ImageFont.truetype(fontloc, 20)
		dtime = astropy.time.TimeDelta(1, format='jd')
		onemin = astropy.time.TimeDelta(60, format='sec')
 
		for div in range(0, (np.size(bminticlocs))): #beta axis minor ticks
			b = d.line([(border+mint,bminticlocs[div]),(border,bminticlocs[div])]
			,fill = (255,255,255,255))
			
		for div in range(0, (np.size(tminticlocs))): #beta axis minor ticks
			b = d.line([(tminticlocs[div],xaxis-mint),(tminticlocs[div],xaxis)]
			,fill = (255,255,255,128))
 
		for div in range(0, (np.size(tmajticlocs))): #simt axis major ticks
			b = d.line([(tmajticlocs[div],xaxis-majt),(tmajticlocs[div],xaxis)],
			fill = (255,255,255,255))
			ticktime = ctime - dtime*tmajticks[div] - onemin
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
		
		if (tlowr <= 0 <= tuppr):
			image_line_x = simt2xpix(0, border, pixwt, tlowr, wscle)
			b = d.line([(image_line_x,border),(image_line_x,2*border+pixhi)]
			,fill = (255,0,0,255))
			d.text((image_line_x - 160,1.25*border), \
			"Time of Image", font=fnt, fill=(255,0,0,255))

		#dustimg.show()
		dustimg.save(image_fixax_img,'png')


