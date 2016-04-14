#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
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
from dust_functions_LINUX import orb_vector, beta2ypix, linsimt2xpix, \
radec_slim

#hard coded parameters - can script these in for softer coding
datafolder = '/unsafe/users/op2/mcnaught/stereo_a/'
fpsav = '/unsafe/users/op2/mcnaught/fpsimsav/'
orbitdir = '/unsafe/users/op2/mcnaught/orbitdata/'
pngsav = '/unsafe/users/op2/mcnaught/dustpngs/'
fitssav = '/unsafe/users/op2/mcnaught/dustfits/'
mapsav = '/unsafe/users/op2/mcnaught/pymapsav/'
txtmaps = '/unsafe/users/op2/mcnaught/txtmappings/'
betau = 2.5
betal = 0.6
bno = 600
simtu = 8.0
simtl = 2.0
tno = 600
comdenom = 'C2006P1'
comname = 'McNaught'
obsloc = 'Stereo A'

#purely for diagnostics
tstart = time.time()

#choosing comet data to us
image_list = sorted(os.listdir(fpsav))
image_total = len(image_list)

#import the orbit data
obsveceq = orb_vector('Stereo_A_orbit_2006_11_21_2007_01_20_xyzvxvyvz_EQ.txt',
                      orbitdir)
comveceq = orb_vector('c2006p1_2006_11_21_2007_01_20_xyzvxvyvz_EQ.txt',
                      orbitdir)
comveceq10 = orb_vector('c2006p1_2005_05_30_2007_01_20_xyzvxvyvz_EQ_10.txt',
                        orbitdir)

for image_id in range(50, 51): #image_total):
    
	image_procname = image_list[image_id].split('.')[0]
	image_basename = image_procname[:21]
	image_simressave = os.path.join(fpsav, image_list[image_id])
	image_fits = os.path.join(datafolder,image_basename + ".fts")
	print (image_basename)

	image_mapped = os.path.join(mapsav,image_procname + "_mapped.npy")
	image_txt_mapping = os.path.join(txtmaps,image_procname + "_map.txt")
	image_dust_img = os.path.join(pngsav,image_procname + "_dustimg.png")
	image_dust_fits = os.path.join(fitssav,image_procname + "_dustfits.fits")

	simres = np.load(image_simressave)

	hdulist = fits.open(image_fits)
	bw_color = (hdulist[0].data)

	fitscoords = image_fits
	onedimg = fits.open(fitscoords)
	w = wcs.WCS(onedimg[0].header, key = 'A')

	#make a 2xN array of all pixel locations
	ya = onedimg[0].data.shape[0]
	xa = onedimg[0].data.shape[1]
	xv, yv = np.meshgrid(np.arange(xa), np.arange(ya))
	xv = np.reshape(xv, (xa*ya,1))
	yv = np.reshape(yv, (xa*ya,1))
	coords = np.zeros((xa*ya,2))
	coords[:,0] = xv[:,0]
	coords[:,1] = yv[:,0]

	#convert each pixel locaion to an RA and DEC array
	radecs = w.wcs_pix2world(coords, 1)
	ra = np.reshape(radecs[:,0], (ya,xa))
	dec = np.reshape(radecs[:,1], (ya,xa))

	#find minimum/maximum values
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

	#get image time from filename
	csec = int(image_basename[13:15])
	cmin = int(image_basename[11:13])
	chour = int(image_basename[9:11])
	cday = int(image_basename[6:8])
	cmonth = int(image_basename[4:6])
	cyear = int(image_basename[0:4])
	ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
				                               chour , cmin, csec))
				                               
	mapping_exists = os.path.exists(image_mapped)
	if (mapping_exists == True): #load if it has
		srcolors = np.load(image_mapped)
	   
	elif (mapping_exists == False): #do if it hasnt
		srcolors = np.zeros((tno,bno,3),dtype=int)
		
		with open(image_txt_mapping, "w") as text_file: 
		
			sync_int_list = np.zeros((tno-1), dtype = int)
			
			for tval in range(0,tno-1):
			
				sync_int_list[tval] = ((int(np.sum(simres[tval,:,12]) >= 1) + 
												int(np.sum(simres[tval+1,:,12]) >= 1))
												== 2)
												
			for ta in np.nonzero(sync_int_list)[0].tolist():
			
				balist1 = np.where(simres[ta+1,:,12] == 1)[0][:-1]
				balist0 = np.where(simres[ta,:,12] == 1)[0][:-1]
				balist = np.intersect1d(balist0,balist1).tolist()
					
				[ra_ta, dec_ta, ra_ta_d0, ra_ta_d1] = radec_slim(ra_m,
									dec, simres[ta,:,10], simres[ta,:,11])
								
				rashape1 = np.shape(ra_ta)[1]
				
				for ba in balist:
						
					text_file.write(('%.3g ,%.3g ,') % 
						(simres[ta,ba,0] , simres[ta,ba,1]))
							
					boxramin = min(simres[ta,ba,10],
										simres[ta+1,ba,10],
				                  simres[ta,ba+1,10],
				                  simres[ta+1,ba+1,10])
				                         	
					boxdemin = min(simres[ta,ba,11],
										simres[ta+1,ba,11],
										simres[ta,ba+1,11],
										simres[ta+1,ba+1,11])
										
					boxramax = max(simres[ta,ba,10],
										simres[ta+1,ba,10],
										simres[ta,ba+1,10],
										simres[ta+1,ba+1,10])
										
					boxdemax = max(simres[ta,ba,11],
										simres[ta+1,ba,11],
                              simres[ta,ba+1,11],
                              simres[ta+1,ba+1,11])
                                   			
					ralocs = np.where((ra_ta > boxramin) &
                    							(ra_ta < boxramax))   
					delocs = np.where((dec_ta > boxdemin) &
                   				 				(dec_ta < boxdemax))
                   				 				
					ralocs1d = ralocs[0]*rashape1+ralocs[1]
					delocs1d = delocs[0]*rashape1+delocs[1]   

					boxlocs1d = np.intersect1d(ralocs1d,delocs1d)
					numin = np.size(boxlocs1d)
							
					if (numin > 0):

						boxlocs = np.zeros((numin,5),dtype = float)
						boxlocs[:,0] = boxlocs1d//rashape1
						boxlocs[:,1] = boxlocs1d%rashape1
								
						boxpath = np.array(
							[[simres[ta,ba,10], simres[ta,ba,11]],
                      [simres[ta+1,ba,10], simres[ta+1,ba,11]],
                      [simres[ta+1,ba+1,10], simres[ta+1,ba+1,11]],
						    [simres[ta,ba+1,10], simres[ta,ba+1,11]],
						    [simres[ta,ba,10], simres[ta,ba,11]]])
                        				
						bbPath = mplPath.Path(boxpath)
                        				
						for n in range(0,numin):
								
							boxlocs[n,2] = ra_ta[int(boxlocs[n,0]),int(boxlocs[n,1])]
							boxlocs[n,3] = dec_ta[int(boxlocs[n,0]),int(boxlocs[n,1])]
							boxlocs[n,4] = bbPath.contains_point((boxlocs[n,2],
                                                                  boxlocs[n,3]))                                        
						boxlocs = boxlocs[np.where(boxlocs[:,4] == 1)]
						numin = np.shape(boxlocs)[0]
							
					if (numin > 0):
					
						text_file.write('Average Enclosed' + ' ,')
						bwtot = 0
						
						for n in range(0,numin):
							bwtot += bw_color[int(boxlocs[n,0])+ra_ta_d0,
													int(boxlocs[n,1])+ra_ta_d1]
													
							text_file.write(str(boxlocs[n,0]+ra_ta_d0) + ' ,' +
                                     str(boxlocs[n,1]+ra_ta_d1) + ' ,')     
                                                       
						srcolors[ta,ba,0] = int(round(float(bwtot)/numin))
						
					else:
						text_file.write('Nearest Neighbour' + ' ,')
						avera = (simres[ta,ba,10] + simres[ta,ba+1,10] +
                           simres[ta+1,ba,10] + simres[ta+1,ba+1,10])/4
						avedec = (simres[ta,ba,11] + simres[ta,ba+1,11] +
                            simres[ta+1,ba,11] + simres[ta+1,ba+1,11])/4       
						distarr = abs(ra_ta - avera) + abs(dec_ta - avedec)
						loc = np.where(distarr == np.min(distarr))
						srcolors[ta,ba,0] = bw_color[loc[0][0]+ra_ta_d0,
																loc[1][0]+ra_ta_d1]
						text_file.write(str(loc[0][0] + ra_ta_d0) + ' ,' +
												str(loc[1][0] + ra_ta_d1) + ' ,')
												
					srcolors[ta,ba,1] = numin
					srcolors[ta,ba,2] = 1
					text_file.write('\n')
               
		np.save(image_mapped,srcolors)

	fits_arr = srcolors[:,:,0].T #indexed by beta first then ejec_t

	if not os.path.exists(image_dust_fits):
		hdu = fits.PrimaryHDU(fits_arr)    
		fitshdr = fits.Header()
		fitshdr['COMMENT'] = "See initialisation file for beta / sim_t values"
		hduhdr = fits.PrimaryHDU(header=fitshdr)
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(image_dust_fits)
  
	if not os.path.exists(image_dust_img):         
	    
		t2sfu = float('%.2g' % simtu)
		t2sfl = float('%.2g' % simtl)
		b1sfu = float('%.1g' % betau)
		b1sfl = float('%.1g' % betal)
    
		pixhi = 1200
		pixwt = int(round(float(pixhi)/bno*tno))
		border = 150
		hscle = pixhi/(np.log10(betau) - np.log10(betal))
		wscle = pixwt/(simtu - simtl)

		dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
				                  	pixhi+int(3*border)),(0,0,0,255))
		d = ImageDraw.Draw(dustimg)

		imgmax = srcolors[:,:,0].max()
		if imgmax < 255: imgmax = 255
        
		fudgefactor = 0.8
		
		sync_int_list = np.zeros((tno-1), dtype = int)
			
		for tval in range(0,tno-1):
			
			sync_int_list[tval] = ((int(np.sum(simres[tval,:,12]) >= 1) + 
												int(np.sum(simres[tval+1,:,12]) >= 1)) == 2)
        
		for ta in np.nonzero(sync_int_list)[0].tolist():
		
			balist1 = np.where(simres[ta+1,:,12] == 1)[0][:-1]
			balist0 = np.where(simres[ta,:,12] == 1)[0][:-1]
			balist = np.intersect1d(balist0,balist1).tolist()
		
			for ba in balist:
			
				fillco = int(round(srcolors[ta,ba,0]*255/imgmax/fudgefactor))
				fillco = sorted([0, fillco, 255])[1]
				b1 = beta2ypix(simres[ta,ba,1], border, pixhi, b1sfl, hscle)
				t1 = linsimt2xpix(simres[ta,ba,0], border, t2sfl, wscle)
				b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, b1sfl, hscle)
				t2 = linsimt2xpix(simres[ta,ba+1,0], border, t2sfl, wscle)
				b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, b1sfl, hscle)
				t3 = linsimt2xpix(simres[ta+1,ba+1,0], border, t2sfl, wscle)
				b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, b1sfl, hscle)
				t4 = linsimt2xpix(simres[ta+1,ba,0], border, t2sfl, wscle)
				a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
				,fill=(fillco,fillco,fillco,255))
             
		a = d.polygon([(border,border),(border*2+pixwt,border), \
        (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
        outline = (255,255,255,128))
    
		bminticks = np.array([0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.6,1.7,
									 1.8,1.9,2.1,2.2,2.3,2.4])
		bmajticks = np.array([0.6,1.0,1.5,2,2.5])
		tmajticks = np.array([2,3,4,5,6,7,8])
		
		bminticlocs = beta2ypix(bminticks, border, pixhi, b1sfl, hscle)
		bmajticlocs = beta2ypix(bmajticks, border, pixhi, b1sfl, hscle)
		tmajticlocs = linsimt2xpix(tmajticks, border, t2sfl, wscle)
 
		majt = 20  #major tick length
		mint = 10  #minor tick length
		xaxis = pixhi + border*2
		fontloc = r'/home/op2/python/lucon.ttf'
		fnt = ImageFont.truetype(fontloc, 30)
		dtime = astropy.time.TimeDelta(1, format='jd')
 
		for div in range(0, (np.size(bminticlocs))): #beta axis minor ticks
			b = d.line([(border+mint,bminticlocs[div]),(border,bminticlocs[div])]
			,fill = (255,255,255,128))
 
		for div in range(0, (np.size(tmajticlocs))): #simt axis major ticks
			b = d.line([(tmajticlocs[div],xaxis-majt),(tmajticlocs[div],xaxis)],
			fill = (255,255,255,128))
			ticktime = ctime - dtime*tmajticks[div]
			tick = ticktime.isot.replace('T','\n')[0:16]
			d.text((tmajticlocs[div] - len(tick)*5,xaxis + 10), \
			tick, font=fnt, fill=(255,255,255,128))
 
		for div in range(0, (np.size(bmajticlocs))): #beta axis major ticks
			b = d.line([(border+majt,bmajticlocs[div]),(border,bmajticlocs[div])]
			,fill = (255,255,255,128))
			tick = str(bmajticks[div])
			d.text((border - len(tick)*5 - 70,bmajticlocs[div] - 10 ), \
			tick, font=fnt, fill=(255,255,255,128))
     
		#axis labels
		d.text((1.5*border + pixwt*0.5 - 150,pixhi + 2.7*border), \
		"Date/Time of Ejection", font=fnt, fill=(255,255,255,128))
		d.text((0.25*border - 10,0.75*border), \
		"Beta", font=fnt, fill=(255,255,255,128))

		#plot title
		plttitle = (comdenom + ' ' + comname + ' from ' + 
		obsloc + ' from date: ' + ctime.isot[0:16].replace('T',' at '))
		tfnt = ImageFont.truetype(fontloc, 40)
		d.text((1.5*border + pixwt*0.5 - len(plttitle)*8 - 290,.35*border), \
		plttitle, font=tfnt, fill=(255,255,255,128))

		dustimg.save(image_dust_img,'png')