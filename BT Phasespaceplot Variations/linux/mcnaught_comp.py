#**********************************************
#Program to visualise two dustplots of comet in
#simt and beta spacewith a fixed axis for time
#**********************************************

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
from dust_functions_LINUX import logbeta2ypix, simt2xpix, round_to_base, RoundToSigFigs

#directories
fpsav_1 = '/unsafe/users/op2/mcnaught_stereo/HI-1/fpsimsav/narrow3/'
fpsav_2 = '/unsafe/users/op2/mcnaught_stereo/HI-2/fpsimsav/20/'
fpsav_3 = '/unsafe/users/op2/mcnaught_soho/C3_Clear/fpsimsav/new/'
mapsav_1 = '/unsafe/users/op2/mcnaught_stereo/HI-1-MGN/pymapsav/narrow3/'
mapsav_2a = '/unsafe/users/op2/mcnaught_stereo/HI-2-diff/pymapsav/20/'
mapsav_2b = '/unsafe/users/op2/mcnaught_stereo/HI-2-MGN/pymapsav/'
mapsav_3 = '/unsafe/users/op2/mcnaught_soho/C3_Clear_MGN/pymapsav/new/'
fixaxsav = '/unsafe/users/op2/mcnaught_comp/'

betau_1 = 3.5
betal_1 = 0.2
bno_1 = 1000
simtu_1 = 16.0
simtl_1 = 2.0
tno_1 = 1000

betau_2 = 2.6
betal_2 = 0.8
bno_2 = 1000
simtu_2 = 20.0
simtl_2 = 6.0
tno_2 = 1000  

betau_3 = 2.0
betal_3 = 0.1
bno_3 = 500
simtu_3 = 5.0
simtl_3 = 1.0
tno_3 = 500

betau = 3.5
betal = 0.2

comdenom = 'C2006P1'
comname = 'McNaught'
obsloc = 'Stereo A'

axlotime = astropy.time.Time(datetime.datetime(2007,1,4,0,0,0))
axhitime = astropy.time.Time(datetime.datetime(2007,1,14,0,0,0))

obslotime = astropy.time.Time(datetime.datetime(2007,1,12,0,0,0))
obshitime = astropy.time.Time(datetime.datetime(2007,1,25,3,0,0))

stime = astropy.time.Time(datetime.datetime(2007,1,20,2,1,0))

hour = astropy.time.TimeDelta(3600, format='sec')
nohours = int(((obshitime - obslotime)/hour).value)

'''
Observation date ranges:
1: Stereo_A: 12th 00:01 --> 18th 22:01
2: Stereo_B: 15th 02:01 --> 25th 02:01
3: Soho:     13th 08:53 --> 15th 18:53
4: TOTAL: 	 12th 00:00 --> 25th 02:00
'''

#choosing comet data to us
image_list_1 = sorted(os.listdir(fpsav_1))
image_total_1 = len(image_list_1)

image_list_2 = sorted(os.listdir(fpsav_2))[1:]
image_total_2 = len(image_list_2)

image_list_3 = sorted(os.listdir(fpsav_3))
image_total_3 = len(image_list_3)

image_times = np.full((nohours,4),-1,dtype=float)
image_ctimes = []

for t in range(0,nohours):
	image_times[t,0] = (obslotime + t*hour).jd
	image_ctimes.append(obslotime + t*hour)

for i in range(0,image_total_1):
	csec = int(image_list_1[i][13:15])
	cmin = int(image_list_1[i][11:13])
	chour = int(image_list_1[i][9:11])
	cday = int(image_list_1[i][6:8])
	cmonth = int(image_list_1[i][4:6])
	cyear = int(image_list_1[i][0:4])
	ctime = astropy.time.Time(datetime.datetime(cyear,cmonth,cday,chour,cmin,csec))
	idx = (np.abs(image_times[:,0]-ctime.jd)).argmin()
	image_times[idx,1] = i
	idx = (np.abs(image_times[:,0]-(ctime+hour).jd)).argmin()
	image_times[idx,1] = i

for i in range(0,image_total_2):
	csec = int(image_list_2[i][13:15])
	cmin = int(image_list_2[i][11:13])
	chour = int(image_list_2[i][9:11])
	cday = int(image_list_2[i][6:8])
	cmonth = int(image_list_2[i][4:6])
	cyear = int(image_list_2[i][0:4])
	ctime = astropy.time.Time(datetime.datetime(cyear,cmonth,cday,chour,cmin,csec))
	idx = (np.abs(image_times[:,0]-ctime.jd)).argmin()
	image_times[idx,2] = i
	idx = (np.abs(image_times[:,0]-(ctime+hour).jd)).argmin()
	image_times[idx,2] = i

	
for i in range(0,image_total_3):
	csec = int(image_list_3[i][13:15])
	cmin = int(image_list_3[i][11:13])
	chour = int(image_list_3[i][9:11])
	cday = int(image_list_3[i][6:8])
	cmonth = int(image_list_3[i][4:6])
	cyear = int(image_list_3[i][0:4])
	ctime = astropy.time.Time(datetime.datetime(cyear,cmonth,cday,chour,cmin,csec))
	idx = (np.abs(image_times[:,0]-ctime.jd)).argmin()
	image_times[idx,3] = i

image_tt = image_times>=0

low1 = -0.7; hih1 = 1.25
low3 = -0.4; hih3 = 0.8
		
pixhi = 600
pixwt = 1000
border = 100

for h in range(0,1):#%START%,%FINISH%):

	ctime = image_ctimes[h]
	imagename = ctime.isot
	image_fixax_img = os.path.join(fixaxsav,imagename + "_fixedcomp.png")
	
	tuppr = ctime.jd - axlotime.jd
	tlowr = ctime.jd - axhitime.jd

	hscle = pixhi/np.log(betau/betal)
	wscle = pixwt/(axhitime.jd - axlotime.jd)

	dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
			                  	pixhi+int(3*border)),(0,0,0,255))
	d = ImageDraw.Draw(dustimg)
	'''
	a = d.polygon([(border,border),(border*2+pixwt,border), \
     (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
     outline = (255,255,255,255))

	tdivmajors = np.array([1.,2.,3.])
	tdivnos = (1/tdivmajors) * (tuppr - tlowr)
	tnodivs = 6
	tdividx = np.where((tdivnos <= tnodivs)==True)[0][0]
	tmajticks = np.arange(tuppr, tlowr-0.0001, -tdivmajors[tdividx])

	tdivminors = np.array([0.5,1,1])
	tminticks = np.arange(tuppr, tlowr-0.0001, -tdivminors[tdividx])
	tminticks = np.setdiff1d(tminticks,tmajticks)
	
	bdivmajors = np.array([0.1,0.2,0.5,1,2])
	bdivnos = (1/bdivmajors) * (betau - betal)
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
	
	for div in range(0, (np.size(bminticlocs))): #beta axis minor ticks
		b = d.line([(border+mint,bminticlocs[div]),(border,bminticlocs[div])]
		,fill = (255,255,255,255))
		
	for div in range(0, (np.size(tminticlocs))): #beta axis minor ticks
		b = d.line([(tminticlocs[div],xaxis-mint),(tminticlocs[div],xaxis)]
		,fill = (255,255,255,128))

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
	
	if (tlowr <= 0 <= tuppr):
		image_line_x = simt2xpix(0, border, pixwt, tlowr, wscle)
		b = d.line([(image_line_x,border),(image_line_x,2*border+pixhi)]
		,fill = (255,0,0,255))
		d.text((image_line_x - 160,1.25*border), \
		"Time of Image", font=fnt, fill=(255,0,0,255))
	'''
	if image_tt[h,2] == True:
	
		image_simres = image_list_2[int(image_times[h,2])]
		image_basename = image_list_2[int(image_times[h,2])].split('.')[0]
		image_srcolor = image_basename + '_mapped.npy'
		
		isec = int(image_basename[13:15])
		imin = int(image_basename[11:13])
		ihour = int(image_basename[9:11])
		iday = int(image_basename[6:8])
		imonth = int(image_basename[4:6])
		iyear = int(image_basename[0:4])
		itime = astropy.time.Time(datetime.datetime(iyear,imonth,iday,ihour,imin,isec))
		
		switchval = (itime.jd - stime.jd)
		image_simres = os.path.join(fpsav_2,image_simres)
		
		if switchval <= 0:
			low2 = -0.3; hih2 = 0.5
			image_srcolor = os.path.join(mapsav_2b,image_srcolor)
			
		if switchval > 0:
			low2 = - 5000 + 900*sorted((0,switchval,999))[1]
			hih2 = 5000 - 900*sorted((0,switchval,999))[1]
			image_srcolor = os.path.join(mapsav_2a,image_srcolor)
		
		simres = np.load(image_simres)
		srcolors = np.load(image_srcolor)

		tcor = 0
		if (itime - ctime).value < -0.02: tcor = 1.0/24

		srcolors_slim = np.copy(srcolors[:,:,0])
		srcolors_slim[np.where((simres[:,:,0] + tcor) >= tuppr)] = 0
		srcolors_slim[np.where((simres[:,:,0] + tcor) <= tlowr)] = 0

		good_bool = srcolors_slim[1:,1:] + srcolors_slim[:-1,:-1] + srcolors_slim[1:,:-1] + srcolors_slim[:-1,1:]
		good_locs = np.where(good_bool==4)
		
		sr_fill = np.clip(np.clip(255.0/(hih2-low2)*(srcolors[:,:,1]-low2),0,255).astype(int),0,255)

		for x in range(0, np.size(good_locs[0])):
			ta = good_locs[0][x]; ba = good_locs[1][x]
			b1 = logbeta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
			t1 = simt2xpix(simres[ta,ba,0] + tcor, border, pixwt, tlowr, wscle)
			b2 = logbeta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
			t2 = simt2xpix(simres[ta,ba+1,0] + tcor, border, pixwt, tlowr, wscle)
			b3 = logbeta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
			t3 = simt2xpix(simres[ta+1,ba+1,0] + tcor, border, pixwt, tlowr, wscle)
			b4 = logbeta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
			t4 = simt2xpix(simres[ta+1,ba,0] + tcor, border, pixwt, tlowr, wscle)
			a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
			,fill=(sr_fill[ta,ba],sr_fill[ta,ba],sr_fill[ta,ba],255))
	
	if image_tt[h,1] == True:
	
		image_simres = image_list_1[int(image_times[h,1])]
		image_basename = image_list_1[int(image_times[h,1])].split('.')[0]
		image_srcolor = image_basename + '_mapped.npy'
		
		image_simres = os.path.join(fpsav_1,image_simres)
		image_srcolor = os.path.join(mapsav_1,image_srcolor)
		
		simres = np.load(image_simres)
		srcolors = np.load(image_srcolor)
			
		isec = int(image_basename[13:15])
		imin = int(image_basename[11:13])
		ihour = int(image_basename[9:11])
		iday = int(image_basename[6:8])
		imonth = int(image_basename[4:6])
		iyear = int(image_basename[0:4])
		itime = astropy.time.Time(datetime.datetime(iyear,imonth,iday,ihour,imin,isec))
		
		tcor = 0
		if (itime - ctime).value < -0.02: tcor = 1.0/24
		
		srcolors_slim = np.copy(srcolors[:,:,0])
		srcolors_slim[np.where((simres[:,:,0] + tcor) >= tuppr)] = 0
		srcolors_slim[np.where((simres[:,:,0] + tcor) <= tlowr)] = 0

		good_bool = srcolors_slim[1:,1:] + srcolors_slim[:-1,:-1] + srcolors_slim[1:,:-1] + srcolors_slim[:-1,1:]
		good_locs = np.where(good_bool==4)
		
		sr_fill = np.clip(255.0/(hih1-low1)*(srcolors[:,:,1]-low1),0,255).astype(int)

		for x in range(0, np.size(good_locs[0])):
			ta = good_locs[0][x]; ba = good_locs[1][x]
			b1 = logbeta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
			t1 = simt2xpix(simres[ta,ba,0] + tcor, border, pixwt, tlowr, wscle)
			b2 = logbeta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
			t2 = simt2xpix(simres[ta,ba+1,0] + tcor, border, pixwt, tlowr, wscle)
			b3 = logbeta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
			t3 = simt2xpix(simres[ta+1,ba+1,0] + tcor, border, pixwt, tlowr, wscle)
			b4 = logbeta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
			t4 = simt2xpix(simres[ta+1,ba,0] + tcor, border, pixwt, tlowr, wscle)
			a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
			,fill=(sr_fill[ta,ba],sr_fill[ta,ba],sr_fill[ta,ba],255))
	
	if image_tt[h,3] == True:
	
		image_simres = image_list_3[int(image_times[h,3])]
		image_basename = image_list_3[int(image_times[h,3])].split('.')[0]
		image_srcolor = image_basename + '_mapped.npy'
		
		image_simres = os.path.join(fpsav_3,image_simres)
		image_srcolor = os.path.join(mapsav_3,image_srcolor)
		
		simres = np.load(image_simres)
		srcolors = np.load(image_srcolor)
		
		srcolors_slim = np.copy(srcolors[:,:,0])
		srcolors_slim[np.where(simres[:,:,0] > tuppr)] = 0
		srcolors_slim[np.where(simres[:,:,0] < tlowr)] = 0
		srcolors_slim[np.where(simres[:,:,1] > betau)] = 0
		srcolors_slim[np.where(simres[:,:,1] < betal)] = 0

		good_bool = srcolors_slim[1:,1:] + srcolors_slim[:-1,:-1] + srcolors_slim[1:,:-1] + srcolors_slim[:-1,1:]
		good_locs = np.where(good_bool==4)
		
		sr_fill = np.clip(255.0/(hih3-low3)*(srcolors[:,:,1]-low3),0,255).astype(int)

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
	
	#dustimg.save(image_fixax_img,'png')
	dustimg.save("base_diff_image_12th_0000",'png')
