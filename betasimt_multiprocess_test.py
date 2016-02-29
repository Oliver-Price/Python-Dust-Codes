#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************
#MULTIPROCESSING TEST
#********************

import string
import easygui
import os
import sys
import pickle
import numpy as np
import matplotlib.path as mplPath
from orbitdata_loading_functions import orb_vector, orb_obs
from plot_functions import beta2ypix, linsimt2xpix, logsimt2xpix, radec_slim, \
greyscale_remap
from astropy.io import fits
import multiprocessing.pool
import time

#%%**********************
#FIRST CELL - GET DATA IN
#************************
#protection against multiprocessing
if __name__ == '__main__': 
    
    #choosing comet data to use
    inputfilefolder = "C:\PhD\Comet_data\Input_files\*pt1.txt"
    inputfile = easygui.fileopenbox(default = inputfilefolder)
    
    #reading main comet parameters
    with open(inputfile, "r") as c:
        cdata = c.readlines()
        comname = cdata[30][12:]
        comdenom = cdata[31][13:-2]
        imagedir = cdata[24][18:-2]
        orbitdir = cdata[25][23:-2]
        idlsav = cdata[26][25:-2]
        pysav = cdata[27][24:-2]
        obsloc = cdata[34][19:24]
        horiztag = cdata[40][10:]
    
    #import the orbit data
    obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                          horiztag, opts = 'obs,eq')
    comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                          horiztag, opts = 'eq')
    comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                            horiztag, opts = 'eq,d10')
    comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)
    
    #choosing fits file to display and getting pathnames
    fitsin = easygui.fileopenbox(default = os.path.join(imagedir, r"*.fits"))
    fitsinfile = os.path.basename(fitsin)
    filebase = fitsinfile[:string.find(fitsinfile,'.')]
    
    #parameter savefile location
    picklesavefile = os.path.join(pysav, filebase + '_dustplot')
    picklexists = os.path.exists(picklesavefile)
    
    greyscale = True
    
    #check vital information exists
    if (picklexists == False):
        sys.exit("Image not calculated")
    
    #import important information
    with open(picklesavefile) as f:
            dparameters = pickle.load(f)
            comcel = dparameters[0]
            comcel10 = dparameters[1]
            rao = dparameters[2]
            deo = dparameters[3]
            rapixl = dparameters[4]
            depixl = dparameters[5]
            ramax = dparameters[6]
            decmax = dparameters[7]
            ramin = dparameters[8]
            decmin = dparameters[9]
            ctime = dparameters[14]
            dtmin = dparameters[15]
            ra = dparameters[16]
            dec = dparameters[17]
            colr = dparameters[18]
            colg = dparameters[19]
            colb = dparameters[20]
            
    #find and import simulation results
    simresdir = os.path.join(pysav, 'simres')
    simin = easygui.fileopenbox(default = os.path.join(simresdir, filebase + '*'))
    simres = np.load(simin)
    
    with open(simin[:-4] + '_parameters') as f:
        sparameters = pickle.load(f)
        tmax = sparameters[0]
        tno = sparameters[2]
        bno = sparameters[3]
        tspace = sparameters[4]

#%%******************************************************
#SECOND CELL - SETTING UP STUFF TO DO THE MULTIPROCESSING
#********************************************************

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def loop(ta, simres = simres, ra = ra, dec = dec, colr = colr, colg = colg, colb = colb):
    [ra_ta, dec_ta, ra_ta_d0, ra_ta_d1] = radec_slim(ra, dec,
                                            simres[ta,:,10], simres[ta,:,11])
    rashape1 = np.shape(ra_ta)[1]
    ba_vals = np.where(simres[ta+1,:,15] == 1)[0][:-1]
    srcolors = np.zeros((np.size(ba_vals),5))
    for ba in ba_vals.tolist():
        boxramin = min(simres[ta,ba,10],simres[ta+1,ba,10],
                       simres[ta,ba+1,10],simres[ta+1,ba+1,10])
        boxdemin = min(simres[ta,ba,11],simres[ta+1,ba,11],
                       simres[ta,ba+1,11],simres[ta+1,ba+1,11])
        boxramax = max(simres[ta,ba,10],simres[ta+1,ba,10],
                       simres[ta,ba+1,10],simres[ta+1,ba+1,10])
        boxdemax = max(simres[ta,ba,11],simres[ta+1,ba,11],
                       simres[ta,ba+1,11],simres[ta+1,ba+1,11])
        ralocs = np.where((ra_ta > boxramin) & (ra_ta < boxramax))   
        delocs = np.where((dec_ta > boxdemin) & (dec_ta < boxdemax))
        ralocs1d = ralocs[0]*rashape1+ralocs[1]
        delocs1d = delocs[0]*rashape1+delocs[1]   
        boxlocs1d = np.intersect1d(ralocs1d,delocs1d)
        numin = np.size(boxlocs1d)
        if (numin > 0): 
            boxlocs = np.empty((numin,5),dtype = float)
            boxlocs[:,0] = np.floor(boxlocs1d/rashape1)
            boxlocs[:,1] = boxlocs1d%rashape1
            boxpath = np.array(
            [[simres[ta,ba,10], simres[ta,ba,11]],
            [simres[ta+1,ba,10], simres[ta+1,ba,11]],
            [simres[ta+1,ba+1,10], simres[ta+1,ba+1,11]],
            [simres[ta,ba+1,10], simres[ta,ba+1,11]]])
            bbPath = mplPath.Path(boxpath)
            for n in xrange(0,numin):
                boxlocs[n,2] = ra_ta[boxlocs[n,0],boxlocs[n,1]]
                boxlocs[n,3] = dec_ta[boxlocs[n,0],boxlocs[n,1]]
                boxlocs[n,4] = bbPath.contains_point((boxlocs[n,2],
                                                      boxlocs[n,3]))                                        
            boxlocs = boxlocs[np.where(boxlocs[:,4] == 1)]
            numin = np.shape(boxlocs)[0]
        if (numin > 0):
            rtot = 0; gtot = 0; btot = 0
            for n in xrange(0,numin):
                rtot += colr[boxlocs[n,0]+ra_ta_d0,boxlocs[n,1]+ra_ta_d1]
                gtot += colg[boxlocs[n,0]+ra_ta_d0,boxlocs[n,1]+ra_ta_d1]
                btot += colb[boxlocs[n,0]+ra_ta_d0,boxlocs[n,1]+ra_ta_d1]                 
            srcolors[ba,0] = int(round(float(rtot)/numin))
            srcolors[ba,1] = int(round(float(gtot)/numin))
            srcolors[ba,2] = int(round(float(btot)/numin))
        else:
            avera = (simres[ta,ba,10] + simres[ta,ba+1,10] +
                    simres[ta+1,ba,10] + simres[ta+1,ba+1,10])/4
            avedec = (simres[ta,ba,11] + simres[ta,ba+1,11] +
                    simres[ta+1,ba,11] + simres[ta+1,ba+1,11])/4       
            distarr = abs(ra_ta - avera) + abs(dec_ta - avedec)
            loc = np.where(distarr == np.min(distarr))
            srcolors[ba,0] = colr[loc[0][0]+ra_ta_d0,loc[1][0]+ra_ta_d1]
            srcolors[ba,1] = colg[loc[0][0]+ra_ta_d0,loc[1][0]+ra_ta_d1]
            srcolors[ba,2] = colb[loc[0][0]+ra_ta_d0,loc[1][0]+ra_ta_d1]
        srcolors[ba,3] = numin
        srcolors[ba,4] = 1
        return srcolors

#%%*********************************************************
#THIRD CELL - DOING THE PIXEL VALUE ASSIGNMENTS
#***********************************************************
#protecting against multiprocessing
if __name__ == '__main__':
    a = time.clock() 
    
    pool = MyPool(4)
    print 'start'
    #results = pool.map(loop, np.arange(len(simres[:,0,0])))
    results = pool.map(loop, np.arange(0,1,1))
    print 'done'
    b = time.clock()

