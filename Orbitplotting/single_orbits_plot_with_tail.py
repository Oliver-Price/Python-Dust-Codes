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
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot Variations")

from particle_sim import *
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

if comdenom == 'c2002v1':
    com_vec = orb_vector_new(comdenom, 'Earth', pysav, orbitdir)
    soho_vec = orb_vector_new(comdenom, 'Soho', pysav, orbitdir, opts = 'obs')
    
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
    
    #comet sim
    cyear = 2007
    cmonth = 1
    cday = 21
    chour = 1

    ctime = Time(datetime.datetime(cyear, cmonth, cday, chour, 0, 0))  
    
    comcel = abs(com_vec[:,0] - ctime.jd).argmin()
    
if comdenom == 'c2011l4':
    
    #start
    syear = 2013
    smonth = 3
    sday = 10
    
    stime = Time(datetime.datetime(syear, smonth, sday, 0, 0, 0))  
    
    scel = abs(earth_vec[:,0] - stime.jd).argmin()
    
    #end
    eyear = 2013
    emonth = 3
    eday = 20
    
    etime = Time(datetime.datetime(eyear, emonth, eday, 0, 0, 0))  
    
    ecel = abs(earth_vec[:,0] - etime.jd).argmin()
    
    #comet sim
    cyear = 2013
    cmonth = 3
    cday = 16
    chour = 0
    cmin = 9

    ctime = Time(datetime.datetime(cyear, cmonth, cday, chour, cmin, 0))  
    
    comcel = abs(com_vec[:,0] - ctime.jd).argmin()

if comdenom == 'c2002v1':
    
    #start
    syear = 2003
    smonth = 2
    sday = 16
    
    stime = Time(datetime.datetime(syear, smonth, sday, 0, 0, 0))  
    
    scel = abs(com_vec[:,0] - stime.jd).argmin()
    
    #end
    eyear = 2003
    emonth = 2
    eday = 21
    
    etime = Time(datetime.datetime(eyear, emonth, eday, 0, 0, 0))  
    
    ecel = abs(com_vec[:,0] - etime.jd).argmin()
    
    #comet sim
    cyear = 2003
    cmonth = 2
    cday = 18
    chour = 8
    cmin = 54

    ctime = Time(datetime.datetime(cyear, cmonth, cday, chour, cmin, 0))  
    
    comcel = abs(com_vec[:,0] - ctime.jd).argmin()
    
#%%

if comdenom == 'c2006p1':
    
    betau = 2.5
    betal = 0.5
    
    tmax = 25
    step = 0.1
    
    tlist = np.arange(step,tmax,step)

    nperday = 50
    output = np.zeros((tlist.size,2,6))
    
    for i, t in np.ndenumerate(tlist):
        cel = abs(com_vec[:,0] - (ctime.jd-t)).argmin()
        start = com_vec[cel,6:12]
        output[i[0],0,:] = part_sim_fine_basic(betal, t, nperday, start)
        output[i[0],1,:] = part_sim_fine_basic(betau, t, nperday, start)
        
    tui = abs(tlist-13).argmin()
    
if comdenom == 'c2011l4':
    
    betau = 3.0
    betal = 0.4
    
    tmax = 14
    step = 0.1
    
    tlist = np.arange(step,tmax,step)

    nperday = 50
    output = np.zeros((tlist.size,2,6))
    
    for i, t in np.ndenumerate(tlist):
        cel = abs(com_vec[:,0] - (ctime.jd-t)).argmin()
        start = com_vec[cel,6:12]
        output[i[0],0,:] = part_sim_fine_basic(betal, t, nperday, start)
        output[i[0],1,:] = part_sim_fine_basic(betau, t, nperday, start)
        
    tui = abs(tlist-7).argmin()

if comdenom == 'c2002v1':
    
    betau = 1.5
    betal = 0.5
    
    tmax = 3.5
    step = 0.1
    
    tlist = np.arange(step,tmax,step)

    nperday = 50
    output = np.zeros((tlist.size,2,6))
    
    for i, t in np.ndenumerate(tlist):
        cel = abs(com_vec[:,0] - (ctime.jd-t)).argmin()
        start = com_vec[cel,6:12]
        output[i[0],0,:] = part_sim_fine_basic(betal, t, nperday, start)
        output[i[0],1,:] = part_sim_fine_basic(betau, t, nperday, start)
        
    tui = abs(tlist-2.5).argmin()

#%%

if comdenom == 'c2006p1':
    
    f0, ax0 = plt.subplots(1,1, figsize=(6,6))
    
    size = 100
    ls = 15
    fs = 20
    d = 0.05
    
    ax0.plot(stereoa_vec[scel:ecel,6],stereoa_vec[scel:ecel,7],c='r',lw=2)
    ax0.scatter(stereoa_vec[ecel,6],stereoa_vec[ecel,7],c='r',marker='x',s=size)
    ax0.plot(stereob_vec[scel:ecel,6],stereob_vec[scel:ecel,7],c='m',lw=2)
    ax0.scatter(stereob_vec[ecel,6],stereob_vec[ecel,7],c='m',marker='x',s=size)
    ax0.plot(soho_vec[scel:ecel,6],soho_vec[scel:ecel,7],c='g',lw=2)
    ax0.scatter(soho_vec[ecel,6],soho_vec[ecel,7],c='g',marker='x',s=size)
    ax0.plot(earth_vec[scel:ecel,6],earth_vec[scel:ecel,7],c='b',lw=2)
    ax0.scatter(earth_vec[ecel,6],earth_vec[ecel,7],c='b',marker='x',s=size)
    ax0.plot(com_vec[scel:ecel,6],com_vec[scel:ecel,7],c='k',lw=2)
    ax0.scatter(com_vec[ecel,6],com_vec[ecel,7],c='k',marker='x',s=size)
    ax0.plot(output[:,0,0],output[:,0,1],'k--')
    ax0.plot(output[:tui,1,0],output[:tui,1,1],'k--')
    ax0.scatter(0,0,c='orange',s=200)
    ax0.text(com_vec[ecel,6]-8.5*d,com_vec[ecel,7]-d,"C/2006 P1",c='k', fontsize=fs)
    ax0.set_xlim([-1.1,0.2])
    ax0.set_ylim([-0.2,1.1])
    ax0.set_title("XY plane (ecliptic)", fontsize = fs)
    ax0.set_xlabel("X (AU)", fontsize = fs)
    ax0.set_ylabel("Y (AU)", fontsize = fs)
    ax0.set_aspect('equal')
    ax0.tick_params(axis='both', which='major', labelsize=ls)
    
    plt.subplots_adjust(left=0.16,
                    bottom=0.11, 
                    right=0.97, 
                    top=0.93,)
    
    plt.savefig('C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\mcnaught_orbits_1.png')
    plt.close()
    
#%%    
    
    f1, ax1 = plt.subplots(1,1, figsize=(6,6))
    
    ax1.plot(stereoa_vec[scel:ecel,7],stereoa_vec[scel:ecel,8],c='r',lw=2)
    ax1.scatter(stereoa_vec[ecel,7],stereoa_vec[ecel,8],c='r',marker='x',s=size)
    ax1.plot(stereob_vec[scel:ecel,7],stereob_vec[scel:ecel,8],c='m',lw=2)
    ax1.scatter(stereob_vec[ecel,7],stereob_vec[ecel,8],c='m',marker='x',s=size)
    ax1.plot(soho_vec[scel:ecel,7],soho_vec[scel:ecel,8],c='g',lw=2)
    ax1.scatter(soho_vec[ecel,7],soho_vec[ecel,8],c='g',marker='x',s=size)
    ax1.plot(earth_vec[scel:ecel,7],earth_vec[scel:ecel,8],c='b',lw=2)
    ax1.scatter(earth_vec[ecel,7],earth_vec[ecel,8],c='b',marker='x',s=size)
    ax1.plot(com_vec[scel:ecel,7],com_vec[scel:ecel,8],c='k',lw=2)
    ax1.scatter(com_vec[ecel,7],com_vec[ecel,8],c='k',marker='x',s=size)
    ax1.plot(output[:,0,1],output[:,0,2],'k--')
    ax1.plot(output[:tui,1,1],output[:tui,1,2],'k--')
    ax1.scatter(0,0,c='orange',s=200)
    ax1.set_xlim([-0.1,1.1])
    ax1.set_ylim([-0.75,0.45])
    ax1.set_title("YZ plane", fontsize = fs)
    ax1.set_xlabel("Y (AU)", fontsize = fs)
    ax1.set_ylabel("Z (AU)", fontsize = fs)
    ax1.set_aspect('equal')
    ax1.tick_params(axis='both', which='major', labelsize=ls)
    
    plt.subplots_adjust(left=0.16,
                    bottom=0.11, 
                    right=0.97, 
                    top=0.93,)
    
    plt.savefig('C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\mcnaught_orbits_2.png')
    plt.close()
    
#%%
    
    f2, ax2 = plt.subplots(1,1, figsize=(6.5,6))
    d = 0.0005
    
    ax2.plot(stereoa_earth_vec[scel:ecel,6],stereoa_earth_vec[scel:ecel,7],c='r',lw=2)
    ax2.scatter(stereoa_earth_vec[ecel,6],stereoa_earth_vec[ecel,7],c='r',marker='x',s=size)
    ax2.plot(stereob_earth_vec[scel:ecel,6],stereob_earth_vec[scel:ecel,7],c='m',lw=2)
    ax2.scatter(stereob_earth_vec[ecel,6],stereob_earth_vec[ecel,7],c='m',marker='x',s=size)
    ax2.plot(soho_earth_vec[scel:ecel,6],soho_earth_vec[scel:ecel,7],c='g',lw=2)
    ax2.scatter(soho_earth_vec[ecel,6],soho_earth_vec[ecel,7],c='g',marker='x',s=size)
    ax2.scatter(0,0,c='b',s=50)
    ax2.text(soho_earth_vec[ecel,6]+d,soho_earth_vec[ecel,7]+d,"SOHO",c='g', fontsize=fs)
    ax2.text(stereob_earth_vec[ecel,6]+2*d,stereob_earth_vec[ecel,7],"STEREO-B",c='m', fontsize=fs)
    ax2.text(stereoa_earth_vec[ecel,6]+d,stereoa_earth_vec[ecel,7]+d,"STEREO-A",c='r', fontsize=fs)
    ax2.text(-0.005+d,-d,"Earth",c='b', fontsize=fs)
    ax2.set_xlim([-0.006,0.018])
    ax2.set_ylim([-0.018,0.006])
    ax2.set_title("XY plane (ecliptic) relative to Earth", fontsize = fs)
    ax2.set_xlabel("X (AU)", fontsize = fs)
    ax2.set_ylabel("Y (AU)", fontsize = fs)
    ax2.set_aspect('equal')
    ax2.tick_params(axis='both', which='major', labelsize=ls)
    
    plt.subplots_adjust(left=0.21,
                    bottom=0.05,
                    right=0.97, 
                    top=1,)
    #plt.subplot_tool()
    
    plt.savefig('C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\mcnaught_orbits_3.png')
    plt.close()
    
#%%

if comdenom == 'c2011l4':
    
    f0, ax0 = plt.subplots(1,1, figsize=(6.5,6))
    size = 100
    ls = 15
    fs = 20
    d = 0.05
    
    ax0.plot(stereoa_vec[scel:ecel,6],stereoa_vec[scel:ecel,7],c='r',lw=2)
    ax0.scatter(stereoa_vec[ecel,6],stereoa_vec[ecel,7],c='r',marker='x',s=size)
    ax0.plot(stereob_vec[scel:ecel,6],stereob_vec[scel:ecel,7],c='m',lw=2)
    ax0.scatter(stereob_vec[ecel,6],stereob_vec[ecel,7],c='m',marker='x',s=size)
    ax0.plot(earth_vec[scel:ecel,6],earth_vec[scel:ecel,7],c='b',lw=2)
    ax0.scatter(earth_vec[ecel,6],earth_vec[ecel,7],c='b',marker='x',s=size)
    ax0.plot(com_vec[scel:ecel,6],com_vec[scel:ecel,7],c='k',lw=2)
    ax0.scatter(com_vec[ecel,6],com_vec[ecel,7],c='k',marker='x',s=size)
    ax0.plot(output[:,0,0],output[:,0,1],'k--')
    ax0.plot(output[:tui,1,0],output[:tui,1,1],'k--')
    ax0.scatter(0,0,c='orange',s=200)
    ax0.text(stereoa_vec[ecel,6]-5*d,stereoa_vec[ecel,7]+2*d,"STEREO-A",c='r', fontsize=fs)
    ax0.text(stereob_vec[ecel,6]-10*d,stereob_vec[ecel,7]+2*d,"STEREO-B",c='m', fontsize=fs)
    ax0.text(earth_vec[ecel,6]+1.5*d,earth_vec[ecel,7]-2*d,"Earth",c='b', fontsize=fs)
    ax0.text(earth_vec[ecel,6]+1.5*d,earth_vec[ecel,7]-2*d,"Earth",c='b', fontsize=fs)    
    ax0.text(com_vec[ecel,6]-13*d,com_vec[ecel,7]+1*d,r"C/2011 L4",c='k', fontsize=fs)    
    ax0.set_xlim([-1.1,1.1])
    ax0.set_ylim([-1.1,1.1])
    ax0.set_title("XY plane (ecliptic)", fontsize = fs)
    ax0.set_xlabel("X (AU)", fontsize = fs)
    ax0.set_ylabel("Y (AU)", fontsize = fs)
    ax0.set_aspect('equal')
    ax0.tick_params(axis='both', which='major', labelsize = ls)

    plt.subplots_adjust(left=0.21,
                    bottom=0.05,
                    right=0.97, 
                    top=1,)
    
    plt.savefig('C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\panstarrs_orbits_1.png')
    plt.close()
    
#%%

    f1, ax1 = plt.subplots(1,1, figsize=(6,6))
      
    ax1.plot(stereoa_vec[scel:ecel,7],stereoa_vec[scel:ecel,8],c='r',lw=2)
    ax1.scatter(stereoa_vec[ecel,7],stereoa_vec[ecel,8],c='r',marker='x',s=size)
    ax1.plot(stereob_vec[scel:ecel,7],stereob_vec[scel:ecel,8],c='m',lw=2)
    ax1.scatter(stereob_vec[ecel,7],stereob_vec[ecel,8],c='m',marker='x',s=size)
    ax1.plot(earth_vec[scel:ecel,7],earth_vec[scel:ecel,8],c='b',lw=2, zorder=1)
    ax1.scatter(earth_vec[ecel,7],earth_vec[ecel,8],c='b',marker='x',s=size)
    ax1.plot(com_vec[scel:ecel,7],com_vec[scel:ecel,8],c='k',lw=2)
    ax1.scatter(com_vec[ecel,7],com_vec[ecel,8],c='k',marker='x',s=size)
    ax1.plot(output[:,0,1],output[:,0,2],'k--')
    ax1.plot(output[:tui,1,1],output[:tui,1,2],'k--')
    ax1.scatter(0,0,c='orange',s=200, zorder=2)
    ax1.set_xlim([-0.9,0.9])
    ax1.set_ylim([-0.9,0.9])
    ax1.set_title("YZ plane", fontsize = fs)
    ax1.set_xlabel("Y (AU)", fontsize = fs)
    ax1.set_ylabel("Z (AU)", fontsize = fs)
    ax1.set_aspect('equal')
    ax1.tick_params(axis='both', which='major', labelsize = ls)

    plt.subplots_adjust(left=0.16,
                    bottom=0.11, 
                    right=0.97, 
                    top=0.94,)
    
    plt.savefig('C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\panstarrs_orbits_2.png')
    plt.close()
     
#%%

if comdenom == 'c2002v1':
    
    f, ax = plt.subplots(1,2, figsize=(13,6))
    size = 100
    ls = 15
    fs = 20
    d = 0.04
    
    ax0.plot(soho_vec[scel:ecel,6],soho_vec[scel:ecel,7],c='g',lw=2)
    ax0.scatter(soho_vec[ecel,6],soho_vec[ecel,7],c='g',marker='x',s=size)
    ax0.plot(com_vec[scel:ecel,6],com_vec[scel:ecel,7],c='k',lw=2)
    ax0.scatter(com_vec[ecel,6],com_vec[ecel,7],c='k',marker='x',s=size)
    ax0.text(soho_vec[ecel,6]+d,soho_vec[ecel,7]-d,"SOHO",c='g', fontsize=fs)
    ax0.text(com_vec[ecel,6]-6*d,com_vec[ecel,7]-2.5*d,r"C/2002 V1",c='k', fontsize=fs)
    ax0.plot(output[:,0,0],output[:,0,1],'k')
    ax0.plot(output[:tui,1,0],output[:tui,1,1],'k')
    ax0.scatter(0,0,c='orange',s=50)
    ax0.set_xlim([-1.0,0.1])
    ax0.set_ylim([-0.3,0.8])
    ax0.set_title("X-Y plane (ecliptic)", fontsize = fs)
    ax0.set_xlabel("X (AU)", fontsize = fs)
    ax0.set_ylabel("Y (AU)", fontsize = fs)
    ax0.set_aspect('equal')
    ax0.tick_params(axis='both', which='major', labelsize = ls)
       
    ax1.plot(soho_vec[scel:ecel,7],soho_vec[scel:ecel,8],c='g',lw=2)
    ax1.scatter(soho_vec[ecel,7],soho_vec[ecel,8],c='g',marker='x',s=size)
    ax1.plot(com_vec[scel:ecel,7],com_vec[scel:ecel,8],c='k',lw=2)
    ax1.scatter(com_vec[ecel,7],com_vec[ecel,8],c='k',marker='x',s=size)
    ax1.plot(output[:,0,1],output[:,0,2],'k--')
    ax1.plot(output[:tui,1,1],output[:tui,1,2],'k--')
    ax1.scatter(0,0,c='orange',s=50)
    ax1.set_xlim([-0.2,0.6])
    ax1.set_ylim([-0.4,0.4])
    ax1.set_title("Y-Z plane", fontsize = fs)
    ax1.set_xlabel("Y (AU)", fontsize = fs)
    ax1.set_ylabel("Z (AU)", fontsize = fs)
    ax1.set_aspect('equal')
    ax1.tick_params(axis='both', which='major', labelsize = ls)
    
    plt.subplots_adjust(left=0.06,
                    bottom=0.105, 
                    right=1.01, 
                    top=0.94, 
                    wspace=0.14, 
                    hspace=0.2)
    
    plt.savefig(r'C:\PhD\Python\Python-Dust-Codes\Orbitplotting\plots\NEAT_orbits_panel.png')