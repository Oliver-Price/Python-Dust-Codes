#*************************************************************
#Program to visualise fits image of comet and overlay a
#finson-probstein diagram according to it's orbital parameters
#*************************************************************

import easygui
import os
import numpy as np
import sys
from astropy.io import fits
from astropy import wcs
from astropy.time import Time
import matplotlib.pyplot as plt

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from orbitdata_loading_functions import *
from io_methods import *

#%%**********************
#FIRST CELL - GET DATA IN
#************************

#choosing comet data to use
inputfilefolder = r"C:\PhD\Comet_data\Input_files\*pt1.txt"
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
    obslocstr = cdata[34][19:]
    horiztag = cdata[40][10:]

#import the orbit data
if comdenom == 'c2006p1':
    com_vec = orb_vector_new(comdenom, 'Earth', pysav, orbitdir)
    soho_vec = orb_vector_new(comdenom, 'Soho', pysav, orbitdir, opts = 'obs')
    stereoa_vec = orb_vector_new(comdenom, 'Stereo_A', pysav, orbitdir, opts = 'obs')
    stereob_vec = orb_vector_new(comdenom, 'Stereo_B', pysav, orbitdir, opts = 'obs')
    earth_vec = orb_vector_new(comdenom, 'Earth', pysav, orbitdir, opts = 'obs')

if comdenom == 'c2011l4':
    com_vec = orb_vector_new(comdenom, 'Earth', pysav, orbitdir)
    stereoa_vec = orb_vector_new(comdenom, 'Stereo_A', pysav, orbitdir, opts = 'obs')
    stereob_vec = orb_vector_new(comdenom, 'Stereo_B', pysav, orbitdir, opts = 'obs')
    earth_vec = orb_vector_new(comdenom, 'Earth', pysav, orbitdir, opts = 'obs')

#%%

if comdenom == 'c2006p1':
    stereoa_earth_vec = stereoa_vec - earth_vec
    stereob_earth_vec = stereob_vec - earth_vec
    soho_earth_vec = soho_vec - earth_vec

#%%

if comdenom == 'c2006p1':
    
    #start
    syear = 2007
    smonth = 1
    sday = 10
    
    stime = Time(datetime.datetime(syear, smonth, sday, 0, 0, 0))  
    
    scel = abs(earth_vec[:,0] - stime.jd).argmin()
    
    #end
    eyear = 2007
    emonth = 1
    eday = 31
    
    etime = Time(datetime.datetime(eyear, emonth, eday, 0, 0, 0))  
    
    ecel = abs(earth_vec[:,0] - etime.jd).argmin()

#%%

if comdenom == 'c2006p1':
    
    f, ax = plt.subplots(1,3, figsize=(20,6))
    size = 50
    
    ax[0].plot(stereoa_vec[scel:ecel,6],stereoa_vec[scel:ecel,7],c='r',lw=2)
    ax[0].scatter(stereoa_vec[ecel,6],stereoa_vec[ecel,7],c='r',marker='x',s=size)
    ax[0].plot(stereob_vec[scel:ecel,6],stereob_vec[scel:ecel,7],c='m',lw=2)
    ax[0].scatter(stereob_vec[ecel,6],stereob_vec[ecel,7],c='m',marker='x',s=size)
    ax[0].plot(soho_vec[scel:ecel,6],soho_vec[scel:ecel,7],c='g',lw=2)
    ax[0].scatter(soho_vec[ecel,6],soho_vec[ecel,7],c='g',marker='x',s=size)
    ax[0].plot(earth_vec[scel:ecel,6],earth_vec[scel:ecel,7],c='b',lw=2)
    ax[0].scatter(earth_vec[ecel,6],earth_vec[ecel,7],c='b',marker='x',s=size)
    ax[0].plot(com_vec[scel:ecel,6],com_vec[scel:ecel,7],c='k',lw=2)
    ax[0].scatter(com_vec[ecel,6],com_vec[ecel,7],c='k',marker='x',s=size)
    ax[0].scatter(0,0,c='orange',s=50)
    ax[0].set_xlim([-1.1,0.1])
    ax[0].set_ylim([-0.1,1.1])
    ax[0].set_title("X-Y plane (ecliptic)", fontsize = 15)
    ax[0].set_xlabel("X (AU)", fontsize = 15)
    ax[0].set_ylabel("Y (AU)", fontsize = 15)
    ax[0].set_aspect('equal')
    ax[0].tick_params(axis='both', which='major', labelsize=13)
    
    ax[1].plot(stereoa_vec[scel:ecel,7],stereoa_vec[scel:ecel,8],c='r',lw=2)
    ax[1].scatter(stereoa_vec[ecel,7],stereoa_vec[ecel,8],c='r',marker='x',s=size)
    ax[1].plot(stereob_vec[scel:ecel,7],stereob_vec[scel:ecel,8],c='m',lw=2)
    ax[1].scatter(stereob_vec[ecel,7],stereob_vec[ecel,8],c='m',marker='x',s=size)
    ax[1].plot(soho_vec[scel:ecel,7],soho_vec[scel:ecel,8],c='g',lw=2)
    ax[1].scatter(soho_vec[ecel,7],soho_vec[ecel,8],c='g',marker='x',s=size)
    ax[1].plot(earth_vec[scel:ecel,7],earth_vec[scel:ecel,8],c='b',lw=2)
    ax[1].scatter(earth_vec[ecel,7],earth_vec[ecel,8],c='b',marker='x',s=size)
    ax[1].plot(com_vec[scel:ecel,7],com_vec[scel:ecel,8],c='k',lw=2)
    ax[1].scatter(com_vec[ecel,7],com_vec[ecel,8],c='k',marker='x',s=size)
    ax[1].scatter(0,0,c='orange',s=50)
    ax[1].set_xlim([-0.1,1.1])
    ax[1].set_ylim([-0.75,0.45])
    ax[1].set_title("Y-Z plane", fontsize = 15)
    ax[1].set_xlabel("Y (AU)", fontsize = 15)
    ax[1].set_ylabel("Z (AU)", fontsize = 15)
    ax[1].set_aspect('equal')
    ax[1].tick_params(axis='both', which='major', labelsize=13)
    
    ax[2].plot(stereoa_earth_vec[scel:ecel,6],stereoa_earth_vec[scel:ecel,7],c='r',lw=2)
    ax[2].scatter(stereoa_earth_vec[ecel,6],stereoa_earth_vec[ecel,7],c='r',marker='x',s=size)
    ax[2].plot(stereob_earth_vec[scel:ecel,6],stereob_earth_vec[scel:ecel,7],c='m',lw=2)
    ax[2].scatter(stereob_earth_vec[ecel,6],stereob_earth_vec[ecel,7],c='m',marker='x',s=size)
    ax[2].plot(soho_earth_vec[scel:ecel,6],soho_earth_vec[scel:ecel,7],c='g',lw=2)
    ax[2].scatter(soho_earth_vec[ecel,6],soho_earth_vec[ecel,7],c='g',marker='x',s=size)
    ax[2].scatter(0,0,c='b',s=50)
    ax[2].set_xlim([-0.006,0.018])
    ax[2].set_ylim([-0.018,0.006])
    ax[2].set_title("X-Y plane (ecliptic) relative to Earth", fontsize = 15)
    ax[2].set_xlabel("X (AU)", fontsize = 15)
    ax[2].set_ylabel("Y (AU)", fontsize = 15)
    ax[2].set_aspect('equal')
    ax[2].tick_params(axis='both', which='major', labelsize=13)
    
    #plt.subplot_tool()
    #plt.show()
    
    plt.subplots_adjust(left=0.035,
                    bottom=0.12, 
                    right=0.99, 
                    top=0.92, 
                    wspace=0.35, 
                    hspace=0.2)

    plt.savefig('C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\mcnaught_orbits.png')