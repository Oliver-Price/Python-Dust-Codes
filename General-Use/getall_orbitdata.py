#*************************************************************
#Program to visualise fits image of comet and overlay a
#finson-probstein diagram according to it's orbital parameters
#*************************************************************

import easygui
import numpy as np
import sys
import astropy.time
import matplotlib.pyplot as plt
import datetime

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from orbitdata_loading_functions import *

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
 
earveceq = orb_vector(comdenom, 'Earth', pysav, orbitdir,
                      horiztag, opts = 'obs')
steraveceq = orb_vector(comdenom, 'Stereo_A', pysav, orbitdir,
                      horiztag, opts = 'obs')
sterbveceq = orb_vector(comdenom, 'Stereo_B', pysav, orbitdir,
                      horiztag, opts = 'obs')
sohoveceq = orb_vector(comdenom, 'Soho', pysav, orbitdir,
                      horiztag, opts = 'obs')

#%%
cyear = 2013
cmonth = 3
cday = 13
chour = 3
cmin = 40

ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                            chour , cmin, 0))  

comcel = np.where(abs(earveceq[:,0] - ctime.jd)==abs(earveceq[:,0] - ctime.jd).min())[0][0]

print(earveceq[comcel,1:6])

a = np.array([earveceq[comcel,6],steraveceq[comcel,6],sterbveceq[comcel,6],sohoveceq[comcel,6]])
b = np.array([earveceq[comcel,7],steraveceq[comcel,7],sterbveceq[comcel,7],sohoveceq[comcel,7]])

plt.scatter(a,b)
plt.axis([-2,2,-2,2])
plt.axes().set_aspect('equal')