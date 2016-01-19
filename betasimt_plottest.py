#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************

import urllib2
import string
import easygui
import os
import sys
import pickle
import numpy as np
import time
import matplotlib.path as mplPath
from orbitdata_loading_functions import orb_vector, orb_obs


#%%

#************************
#FIRST CELL - GET DATA IN
#************************

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

#check if image exists already
imgsav = os.path.join(imagedir, filebase + '.png')
imgexists = os.path.exists(imgsav)

#parameter savefile location
picklesavefile = os.path.join(pysav, filebase + '_dustplot')
picklexists = os.path.exists(picklesavefile)

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
        border = dparameters[10]
        pixheight = dparameters[11]
        scale = dparameters[12]
        ctime = dparameters[13]
        dtmin = dparameters[14]
        ra = dparameters[15]
        dec = dparameters[16]
        colr = dparameters[17]
        colg = dparameters[18]
        colb = dparameters[19]
        
#find and import simulation results
simresdir = os.path.join(pysav, 'simres')
simin = easygui.fileopenbox(default = os.path.join(simresdir, filebase + '*'))
with open(simin) as f:
        sparameters = pickle.load(f)
        simres = sparameters[0]
        tmax = sparameters[1]
        bmax = sparameters[2]
        tno = sparameters[3]
        bno = sparameters[4]
    
#%%

#****************************************************
#SECOND CELL - FIND RELEVANT DATA FROM ORIGINAL IMAGE
#****************************************************
srcolors = np.empty((tno,bno,3),dtype = int)
for ta in xrange(0, tno-1):
    for ba in xrange(0, bmax[ta+1]-1):
        boxramin = min(simres[ta,ba,10],simres[ta+1,ba,10],
                       simres[ta,ba+1,10],simres[ta+1,ba+1,10])
        boxdemin = min(simres[ta,ba,11],simres[ta+1,ba,11],
                       simres[ta,ba+1,11],simres[ta+1,ba+1,11])
        boxramax = max(simres[ta,ba,10],simres[ta+1,ba,10],
                       simres[ta,ba+1,10],simres[ta+1,ba+1,10])
        boxdemax = max(simres[ta,ba,11],simres[ta+1,ba,11],
                       simres[ta,ba+1,11],simres[ta+1,ba+1,11])  
        ralocs = np.where((ra > boxramin) & (ra < boxramax))   
        delocs = np.where((dec > boxdemin) & (dec < boxdemax))
        rashape0 = np.shape(ra)[0]
        ralocs1d = ralocs[0]*rashape0+ralocs[1]
        delocs1d = delocs[0]*rashape0+delocs[1]   
        boxlocs1d = np.intersect1d(ralocs1d,delocs1d)
        numin = np.size(boxlocs1d)
        if (numin > 0): 
            boxlocs = np.empty((numin,3),dtype = int)
            boxlocs[:,0] = np.floor(boxlocs1d/rashape0)
            boxlocs[:,1] = boxlocs1d%rashape0
            bbPath = mplPath.Path(np.array(
            [[simres[ta,ba,10], simres[ta,ba,11]],
            [simres[ta+1,ba,10], simres[ta+1,ba,11]],
            [simres[ta+1,ba+1,10], simres[ta+1,ba+1,11]],
            [simres[ta,ba+1,10], simres[ta,ba+1,11]]]))#
            for n in xrange(0,numin):
                boxlocs[n,2] = bbPath.contains_point((boxlocs[n,0],
                                                        boxlocs[n,1]))
            boxlocs = boxlocs[np.where(boxlocs[:,2] == 1)]
            numin = np.shape(boxlocs)[0]
            if (numin == 1):
                srcolors[ta,ba,0] = colr[boxlocs[0,0],boxlocs[0,1]]
                srcolors[ta,ba,1] = colg[boxlocs[0,0],boxlocs[0,1]]
                srcolors[ta,ba,2] = colb[boxlocs[0,0],boxlocs[0,1]]
            if (numin > 1):
                for n in xrange(0,numin):
                     rtot += colr[boxlocs[n,0],boxlocs[n,1]]
                     gtot += colg[boxlocs[n,0],boxlocs[n,1]]
                     btot += colb[boxlocs[n,0],boxlocs[n,1]]
                srcolors[ta,ba,0] = rtot/numin
                srcolors[ta,ba,1] = gtot/numin
                srcolors[ta,ba,2] = btot/numin   
        if (numin == 0):
            avera = (simres[ta,ba,10] + simres[ta,ba+1,10] +
                    simres[ta+1,ba,10] + simres[ta+1,ba+1,10])/4
            avedec = (simres[ta,ba,11] + simres[ta,ba+1,11] +
                    simres[ta+1,ba,11] + simres[ta+1,ba+1,11])/4       
            distarr = abs(ra - avera) + abs(dec - avedec)
            loc = np.where(distarr == np.min(distarr))
            srcolors[ta,ba,0] = colr[loc[0][0],loc[1][0]]
            srcolors[ta,ba,1] = colg[loc[0][0],loc[1][0]]
            srcolors[ta,ba,2] = colb[loc[0][0],loc[1][0]]
            
            

