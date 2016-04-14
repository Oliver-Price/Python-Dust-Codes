#*************************************************************
#Program to visualise fits image of comet and overlay a
#finson-probstein diagram according to it's orbital parameters
#*************************************************************

import os
import sys
import time
import string
import pickle
import datetime
import numpy as np
from astropy.io import fits
from astropy import wcs
import astropy.time
import matplotlib.path as mplPath
from dust_functions_LINUX import orb_vector, pos2radec, part_sim

#hard coded parameters - can script these in for softer coding
datafolder = '/unsafe/users/op2/mcnaught/stereo_a/'
pysav = '/unsafe/users/op2/mcnaught/fpsimsav/'
orbitdir = '/unsafe/users/op2/mcnaught/orbitdata/'
betau = 2.5
betal = 0.6
bno = 100
simtu = 8.0
simtl = 2.0
tno = 100

tstart = time.time()

#choosing comet data to us
image_list = sorted(os.listdir(datafolder))
image_total = len(image_list)

#import the orbit data
obsveceq = orb_vector('Stereo_A_orbit_2006_11_21_2007_01_20_xyzvxvyvz_EQ.txt',
                      orbitdir)
comveceq = orb_vector('c2006p1_2006_11_21_2007_01_20_xyzvxvyvz_EQ.txt',
                      orbitdir)
comveceq10 = orb_vector('c2006p1_2005_05_30_2007_01_20_xyzvxvyvz_EQ_10.txt',
                        orbitdir)

for image_id in range(75, 76): #image_total):

	image_basename = image_list[image_id].split('.')[0]
	image_fits = os.path.join(datafolder,image_list[image_id])
	print (image_basename)

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

	#find relevant observer parameters of comet at observer time
	comcel = np.where(abs(obsveceq[:,0] - ctime.jd) < 1e-4)[0][0]

	#find closest matching cell in dt = 10min data
	comcel10 = np.where(abs(comveceq10[:,0] - ctime.jd) < 3.5e-3)[0][0]
	c10t = astropy.time.Time(comveceq10[comcel10,0], format = 'jd')

	#find distance to actual time
	dt = ctime - c10t
	dtmin = int(round(dt.sec/60))

	#find ra and dec of comet
	com_ra_dec = pos2radec(comveceq[comcel,6:9] - obsveceq[comcel,6:9])

	#prep for finding image case
	tvals = np.linspace(simtl, simtu, tno)
	bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))
	efinp = obsveceq[comcel,6:9]

	#prepare for raycasting
	com_box_path = np.array(
                    	[[ra_m[0,0],		dec[0,0]		],
                   		 [ra_m[ya-1,0],	 	dec[ya-1,0]		],
                  		 [ra_m[ya-1,xa-1],  dec[ya-1,xa-1]	],
                   		 [ra_m[0,xa-1],	 	dec[0,xa-1]		],
                   		 [ra_m[0,0],		dec[0,0]		]]) 
	com_Path = mplPath.Path(com_box_path)
	
	#check if comet is in image
	com_in_image = com_Path.contains_point((com_ra_dec[0],com_ra_dec[1]))
	if com_in_image == True:
		image_case = 1;
		
	#if it isn't, check for synchrone intersection with image
	else:
		
		sync_intersects = np.zeros(tno,dtype = int)
		
		for tval_test in range(0, tno):	
			simt10min = int(round(144*tvals[tval_test]))
			pstart = comveceq10[comcel10-simt10min,6:12]
			sim_lo = part_sim(bvals[0],simt10min,30,3,pstart,efinp,dtmin)
			sim_hi = part_sim(bvals[bno-1],simt10min,30,3,pstart,efinp,dtmin)
			ra_dec_lo_test = pos2radec(sim_lo[1][0:3] -
						obsveceq[int(comcel-10*simt10min+sim_lo[0]),6:9])
			ra_dec_hi_test = pos2radec(sim_hi[1][0:3] -
						obsveceq[int(comcel-10*simt10min+sim_hi[0]),6:9])
			sync_box_path = np.array(
                    	[[ra_dec_lo_test[0],ra_dec_lo_test[1]],
                   		 [ra_dec_hi_test[0],ra_dec_hi_test[1]]])
			sync_Path = mplPath.Path(sync_box_path)
			sync_intersects[tval_test] = com_Path.intersects_path(sync_Path)
            
		if (np.sum(sync_intersects) > 0):
			image_case = 2
		elif (np.sum(sync_intersects) == 0): image_case = 3
	
	if image_case == 1: #IF COMET IS INSIDE IMAGE
					  
		#prep for simulation loop proper
		simres = np.empty((tno,bno,13),dtype = float)
		tidx = 0
		
		#this is the loop that does the business
		while (tidx < tno):
			bidx = 0
			point_in_image = 1
			while ( bidx < bno and point_in_image == 1):
				simt10min = int(round(144*tvals[tidx]))
				pstart = comveceq10[comcel10-simt10min,6:12]
				sim = part_sim(bvals[bidx],simt10min,30,3,pstart,efinp,dtmin)
				simres[tidx,bidx,0] = float(simt10min)/144
				simres[tidx,bidx,1] = bvals[bidx]
				simres[tidx,bidx,2] = sim[0] #length of simulation in minutes
				simres[tidx,bidx,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
				simres[tidx,bidx,4:10] = sim[1] #finishing pos/vel
				simres[tidx,bidx,10:12] = pos2radec(sim[1][0:3] - 
				obsveceq[int(simres[tidx,bidx,3]),6:9])
				point_in_image = int(com_Path.contains_point(
								 (simres[tidx,bidx,10],simres[tidx,bidx,11])))
				simres[tidx,bidx,12] = point_in_image
				bidx += 1
			tidx += 1
		

		simressavefile = os.path.join(pysav, image_basename + '_' + str(betal) + '_'
				         + str(betau) + '_' + str(bno)+ '_' + str(simtl) + '_'
				         + str(simtu) + '_' + str(tno))
		simressavefile = simressavefile.replace('.','\'')
		np.save(simressavefile, simres)
		
	if image_case == 2: #IF COMET OUTSIDE IMAGE
		
		#prep for simulation loop proper
		simres = np.zeros((tno,bno,13),dtype = float)
		tidx_list = np.where(sync_intersects == 1)[0]
		syn_exit_tt = np.zeros((2,2),dtype = int)
		syn_exit_tt[0,1] = 1
		tidx = 0
		
		#this is the loop that does the business
		for tidx_list_val in range(0,np.size(tidx_list)):
			bidx = bno-1
			tidx = tidx_list[tidx_list_val]
			syn_exited_image = 0
			prev_stat = 0
			while (bidx >= 0 and syn_exited_image == 0):
				simt10min = int(round(144*tvals[tidx]))
				pstart = comveceq10[comcel10-simt10min,6:12]
				sim = part_sim(bvals[bidx],simt10min,30,3,pstart,efinp,dtmin)
				simres[tidx,bidx,0] = float(simt10min)/144
				simres[tidx,bidx,1] = bvals[bidx]
				simres[tidx,bidx,2] = sim[0] #length of simulation in minutes
				simres[tidx,bidx,3] = comcel-10*simt10min+sim[0] #find relevant cell for end of traj
				simres[tidx,bidx,4:10] = sim[1] #finishing pos/vel
				simres[tidx,bidx,10:12] = pos2radec(sim[1][0:3] - 
				obsveceq[int(simres[tidx,bidx,3]),6:9])
				point_in_image = int(com_Path.contains_point(
								 (simres[tidx,bidx,10],simres[tidx,bidx,11])))
				simres[tidx,bidx,12] = point_in_image
				syn_exited_image = syn_exit_tt[point_in_image,prev_stat]
				prev_stat = point_in_image
				bidx -= 1

		simressavefile = os.path.join(pysav, image_basename + '_' + str(betal) + '_'
				         + str(betau) + '_' + str(bno)+ '_' + str(simtl) + '_'
				         + str(simtu) + '_' + str(tno))
		simressavefile = simressavefile.replace('.','\'')
		np.save(simressavefile, simres)
		
	else:
		print("Image " + image_basename + " was no good.")
		
print (time.time() - tstart)