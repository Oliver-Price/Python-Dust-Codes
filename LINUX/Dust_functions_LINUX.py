#LINUX various functions for simulation code

import os
import numpy as np
import math as m

#%% - orbit data loading function

def orb_vector(dataname , datapath):

	savename = dataname.split('.')[0] + '.npy'
	savefile = os.path.join(datapath, savename)
	saveexists = os.path.isfile(savefile)

	if saveexists == False: #if it doesn't exist read data from file
		dataread = os.path.join(datapath, dataname)
		with open(dataread, "r") as dat:
                
			#reads raw data as string 
			rawdata = dat.readlines()
			datasize = int(len(rawdata)/3)
			data = np.empty((datasize,13),dtype = float)
         
			#converts string to numpy array
			for hrow in range(0,datasize,1):
				data[hrow,0] = rawdata[3*hrow][0:17]            #jd
				data[hrow,1] = rawdata[3*hrow][25:29]           #year
				data[hrow,2] = mon2num(rawdata[3*hrow][30:33])  #month
				data[hrow,3] = rawdata[3*hrow][34:36]           #day
				data[hrow,4] = rawdata[3*hrow][37:39]           #hour
				data[hrow,5] = rawdata[3*hrow][40:42]           #minute
				p = np.fromstring(rawdata[(hrow*3)+1],
								 dtype=float, count=3, sep='  ')
				v = np.fromstring(rawdata[(hrow*3)+2],
								 dtype=float, count=3, sep='  ')
				data[hrow,6:9] = p          #positions xyz
				data[hrow,9:12] = v         #velocities vxvyvz
				data[hrow,12] = np.linalg.norm(data[hrow,6:9])
			data.dump(savefile)   #saves data for future loading
            
	elif saveexists == True: #if save exists, load that
		data = np.load(savefile)
		print ('Loading saved data')
	return data
    
#%% - orbit data loading function

def orb_obs(dataname , datapath):

	savename = dataname.split('.')[0] + '.npy'
	savefile = os.path.join(datapath, savename)
	saveexists = os.path.isfile(savefile)
    
	if saveexists == False: #if it doesn't exist read data from file
		dataread = os.path.join(datapath, dataname)
		with open(dataread, "r") as dat:
        
			#reads raw data as string 
			rawdata = dat.readlines()
			data = np.empty((len(rawdata),7),dtype = float)
			
			for row in range(0,len(rawdata),1):
				data[row,0] = rawdata[row][1:5]             #year
				data[row,1] = mon2num(rawdata[row][6:9])    #month
				data[row,2] = rawdata[row][10:12]           #day
				data[row,3] = rawdata[row][13:15]           #hour
				data[row,4] = rawdata[row][16:18]           #minute
				data[row,5] = rawdata[row][23:32]           #ra
				data[row,6] = rawdata[row][33:42]           #dec
			
			data.dump(savefile)
			
	elif saveexists == True: #if save exists, load that
		data = np.load(savefile)
		print ('Loading saved data')
	return data      

#%% - converts month text to a number

def mon2num(month):
    if month == 'Jan':
        mnum = 1
    elif month == 'Feb':
        mnum = 2
    elif month == 'Mar':
        mnum = 3
    elif month == 'Apr':
        mnum = 4
    elif month == 'May':
        mnum = 5
    elif month == 'Jun':
        mnum = 6
    elif month == 'Jul':
        mnum = 7
    elif month == 'Aug':
        mnum = 8
    elif month == 'Sep':
        mnum = 9
    elif month == 'Oct':
        mnum = 10
    elif month == 'Nov':
        mnum = 11
    elif month == 'Dec':
        mnum = 12
    return mnum
    
#%% - converts earth centric equatorial xyz coord to ra and dec
    
def pos2radec(position):
    ra = m.atan(position[1]/position[0])*360/(2*m.pi)
    if position[0] >= 0:
        if position[1] < 0:
            ra = ra + 360
    if position[0] < 0:
        ra = ra + 180
    r = np.linalg.norm(position)
    dec = m.asin(position[2]/r)*360/(2*m.pi)
    return (ra,dec)
   
#%% - RK4 METHOD for simulating dust particle trajectories

'''
INPUTS: BETA, TIME TO SIMULATE FOR, DT,
        INPUT POSITION AND VELOCITIES, POSITION OF EARTH AT FINISHING TIME

simt in 10s of minutes , ndt is now number of dt
ltdt is the dt when calculated lt correction

Diff function does this:
Input = Particle[r,v]	(Array with position and velocity)
dr = v;
a = -GM/|r|*|r|;
dv = a dot R(comet->sun) (Project a radially)
return dParticle = Array[dr,dv];

Distance is in AU with speed in AU/day but the simulation runs in minutes.

For the RK4 method we convert minutes to seconds then:
1.1574074074074073e-05 converts AU/day to AU/seconds

We also need solar GM in units of AU^3/day/second
Solar GM = 1.32712440018e20 m^3 s^-2
3.4249098175024633e-09 = 1.32712440018e20 / 149597870700^3 * 86400

For the output we assume the earth doesn't move a large distance away from
comet at end of it's traj. The final position is given where the simt is
reached accounting for a LT time to earth.
8.316746397269274 is c in Au/minute.

OUTPUTS: Final location and finishing time
'''

def part_sim(beta, simt, ndt, ltdt, pstart, efinp, cor): 
    simt = simt*10 + cor #simt in minutes
    t = 0
    dt = (simt-30)/ndt
    traj = pstart
    for n in range(0,ndt,1): #go to up to 30 minutes before current
        traj = rk4(traj, beta, dt) #traj updated
        t += dt    
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)
    
def diff(pstate, beta):
    diffed = np.empty((6),dtype = float)
    diffed[0:3] = 1.1574074074074073e-05*pstate[3:6]
    invr = 1/(np.linalg.norm(pstate[0:3]))
    acc = (-3.4249098175024633e-09)*(1-beta)*invr*invr*invr
    diffed[3:6] = acc*pstate[0:3] #3rd invr normalises the position vector
    return diffed
    
def rk4(pstate, beta, dt): 
    dt = dt*60 #DT INPUT IN MINUTES
    k1 = diff(pstate, beta)*dt
    k2 = diff(pstate+k1*0.5, beta)*dt
    k3 = diff(pstate+k2*0.5, beta)*dt
    k4 = diff(pstate+k3, beta)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666
    
#%% - Plot methods for making axis for dust plot images
    
def simt2xpix(simt, border, pixwidth, simtl, scale):
    return pixwidth + border*1.5 + (simtl - simt)*scale
    
def beta2ypix(beta, border, pixheight, betal, scale):
    return pixheight + border*1.5 + (betal - beta)*scale

#%% FUNCTION TO SLIM DOWN AN ARRAY OF RA AND DEC VALUES USING ANOTHER SET OF RA
#AND DEC VALUES AS REFRENCE

def radec_slim(ra_data, dec_data, ra_ref, dec_ref):

    #this assumes that the data arrays of ra and dec have the same shape
    rashape0 = np.shape(ra_data)[0]
    rashape1 = np.shape(ra_data)[1]
    
    raimin = np.min(ra_ref[np.nonzero(ra_ref)])
    raimax = np.max(ra_ref[np.nonzero(ra_ref)])
    decimin = np.min(dec_ref[np.nonzero(dec_ref)])
    decimax = np.max(dec_ref[np.nonzero(dec_ref)])
    
    distmin = abs(ra_data - raimin) + abs(dec_data - decimin)
    minloc = np.where(distmin == np.min(distmin))
    distmax = abs(ra_data - raimax) + abs(dec_data - decimax)
    maxloc = np.where(distmax == np.min(distmax))
    
    up_index_0 = max(minloc[0][0], maxloc[0][0])
    down_index_0 = min(minloc[0][0], maxloc[0][0])
    up_index_1 = max(minloc[1][0], maxloc[1][0])
    down_index_1 = min(minloc[1][0], maxloc[1][0])
    
    if up_index_0 != (rashape0 - 1):
        up_index_0 += 1
    if down_index_0 != 0:
        down_index_0 -= 1
    if up_index_1 != (rashape1 - 1):
        up_index_1 += 1
    if down_index_1 != 0:
        down_index_1 -= 1
        
    ra_slim = ra_data[down_index_0:up_index_0,down_index_1:up_index_1]
    dec_slim = dec_data[down_index_0:up_index_0,down_index_1:up_index_1]

    return(ra_slim, dec_slim, down_index_0, down_index_1)
