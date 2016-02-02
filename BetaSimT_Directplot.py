#*************************************************************
#Program to visualise dustplot of comet in simt and beta space
#*************************************************************

import string
import easygui
import os
import sys
import pickle
import numpy as np
import time
import matplotlib.path as mplPath
from orbitdata_loading_functions import orb_vector, orb_obs
from plot_functions import beta2ypix, linsimt2xpix, logsimt2xpix, radec_slim
from PIL import Image, ImageDraw, ImageFont


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
    bmax = sparameters[1]
    tno = sparameters[2]
    bno = sparameters[3]
    tspace = sparameters[4]

#%%

#***********************************************************
#SECOND CELL - IDENTIFYING PIXEL VALUES TO USE FOR BETA/SIMT
#***********************************************************
#fig = plt.figure()
srcolors = np.zeros((tno,bno,5),dtype=int)

for ta in xrange(0, tno-1):
    [ra_ta, dec_ta, ra_ta_d0, ra_ta_d1] = radec_slim(ra, dec,
                                            simres[ta,:,10], simres[ta,:,11])
    rashape1 = np.shape(ra_ta)[1]
    for ba in xrange(0, bmax[ta+1]):
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
            srcolors[ta,ba,0] = rtot/numin
            srcolors[ta,ba,1] = gtot/numin
            srcolors[ta,ba,2] = btot/numin
        else:
            avera = (simres[ta,ba,10] + simres[ta,ba+1,10] +
                    simres[ta+1,ba,10] + simres[ta+1,ba+1,10])/4
            avedec = (simres[ta,ba,11] + simres[ta,ba+1,11] +
                    simres[ta+1,ba,11] + simres[ta+1,ba+1,11])/4       
            distarr = abs(ra - avera) + abs(dec - avedec)
            loc = np.where(distarr == np.min(distarr))
            srcolors[ta,ba,0] = colr[loc[0][0],loc[1][0]]
            srcolors[ta,ba,1] = colg[loc[0][0],loc[1][0]]
            srcolors[ta,ba,2] = colb[loc[0][0],loc[1][0]]
        srcolors[ta,ba,3] = numin
        srcolors[ta,ba,4] = 1
    print float(ta*100)/tno

#%%

#*********************************
#THIRD CELL - PLOT DATA ONTO IMAGE
#*********************************          
            
simtl = simres[0,0,0]; simtu = simres[tno-1,0,0]
betal = simres[0,0,1]; betau = simres[0,bno-1,1]

t1sfu = float('%.1g' % simtu)
t1sfl = float('%.1g' % simtl)
b1sfu = float('%.1g' % betau)
b1sfl = float('%.1g' % betal)

pixhi = 800
pixwt = 1600
border = 100
hscle = pixhi/(np.log10(betau) - np.log10(betal))
if (tspace == 'log'):
    wscle = pixwt/(np.log10(simtu) - np.log10(simtl))
elif (tspace == 'lin'):
    wscle = pixwt/(simtu - simtl)
    
dustimg = Image.new('RGBA', (pixwt+int(2.5*border),
                             pixhi+int(3*border)),(0,0,0,255))
d = ImageDraw.Draw(dustimg)
greyscale = False
nmmax = np.log(np.max(srcolors[:,:,3])+1)

if (greyscale == True):  
    for ta in xrange(0, tno-1):
        for ba in xrange(0, bmax[ta+1]):
            fillco = int(round(0.333333333333333333*(srcolors[ta,ba,0] + 
            srcolors[ta,ba,1] + srcolors[ta,ba,2])))
            b1 = beta2ypix(simres[ta,ba,1], border, pixhi, b1sfl, hscle)
            t1 = linsimt2xpix(simres[ta,ba,0], border, t1sfl, wscle)
            b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, b1sfl, hscle)
            t2 = linsimt2xpix(simres[ta,ba+1,0], border, t1sfl, wscle)
            b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, b1sfl, hscle)
            t3 = linsimt2xpix(simres[ta+1,ba+1,0], border, t1sfl, wscle)
            b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, b1sfl, hscle)
            t4 = linsimt2xpix(simres[ta+1,ba,0], border, t1sfl, wscle)
            a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
            ,fill=(fillco,fillco,fillco,255))
else:
    for ta in xrange(0, tno-1):
        for ba in xrange(0, bmax[ta+1]):
            b1 = beta2ypix(simres[ta,ba,1], border, pixhi, b1sfl, hscle)
            t1 = linsimt2xpix(simres[ta,ba,0], border, t1sfl, wscle)
            b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, b1sfl, hscle)
            t2 = linsimt2xpix(simres[ta,ba+1,0], border, t1sfl, wscle)
            b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, b1sfl, hscle)
            t3 = linsimt2xpix(simres[ta+1,ba+1,0], border, t1sfl, wscle)
            b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, b1sfl, hscle)
            t4 = linsimt2xpix(simres[ta+1,ba,0], border, t1sfl, wscle)
            a = d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)]
            ,fill=(srcolors[ta,ba,0],srcolors[ta,ba,1],srcolors[ta,ba,2],255))
            
#%%

#***********************
#FOURTH CELL - DRAW AXIS
#***********************

a = d.polygon([(border,border),(border*2+pixwt,border), \
    (border*2+pixwt,border*2+pixhi),(border,border*2+pixhi)], \
    outline = (255,255,255,128))

decades = np.logspace(-4,4,9)
tdecl = np.searchsorted(decades,simtl, side = 'right')-1
tdecu = np.searchsorted(decades,simtu)
bdecl = np.searchsorted(decades,betal, side = 'right')-1
bdecu = np.searchsorted(decades,betau)

bminticks = np.linspace(decades[bdecl],decades[bdecl+1],10)
for bdec in xrange(bdecl+1, bdecu):
    bminticks = np.concatenate((bminticks,
                      np.linspace(decades[bdec],decades[bdec+1],10)[1:10]))
bminticks = bminticks[np.searchsorted(bminticks,b1sfl):
                          np.searchsorted(bminticks,b1sfu)+1]
bmajticks = np.intersect1d(bminticks,decades)

tminticks = np.linspace(decades[tdecl],decades[tdecl+1],10)
for tdec in xrange(tdecl+1, tdecu):
    tminticks = np.concatenate((tminticks,
                      np.linspace(decades[tdec],decades[tdec+1],10)[1:10]))
tminticks = tminticks[np.searchsorted(tminticks,t1sfl):
                          np.searchsorted(tminticks,t1sfu)+1]
tmajticks = np.intersect1d(tminticks,decades)

bminticlocs = beta2ypix(bminticks, border, pixhi, b1sfl, hscle)
tminticlocs = linsimt2xpix(tminticks, border, t1sfl, wscle)
bmajticlocs = beta2ypix(bmajticks, border, pixhi, b1sfl, hscle)
tmajticlocs = linsimt2xpix(tmajticks, border, t1sfl, wscle)

majt = 20  #major tick length
mint = 10  #minor tick length
xaxis = pixhi + border*2
fontloc = r'C:\Windows\winsxs\amd64_microsoft-windows-f..etype-lucidaconsole_31bf3856ad364e35_6.1.7600.16385_none_5b3be3e0926bd543\lucon.ttf'
fnt = ImageFont.truetype(fontloc, 20)

for div in xrange(0, (np.size(tminticlocs))): #simt axis minor ticks
    b = d.line([(tminticlocs[div],xaxis-mint),(tminticlocs[div],xaxis)],\
    fill = (255,255,255,128))

for div in xrange(0, (np.size(bminticlocs))): #beta axis minor ticks
    b = d.line([(border+mint,bminticlocs[div]),(border,bminticlocs[div])],\
    fill = (255,255,255,128))

for div in xrange(0, (np.size(tmajticlocs))): #simt axis major ticks
    b = d.line([(tmajticlocs[div],xaxis-majt),(tmajticlocs[div],xaxis)],\
    fill = (255,255,255,128))
    tick = str(tmajticks[div])
    d.text((tmajticlocs[div] - len(tick)*5,xaxis + 10), \
    tick, font=fnt, fill=(255,255,255,128))

for div in xrange(0, (np.size(bmajticlocs))): #beta axis major ticks
    b = d.line([(border+majt,bmajticlocs[div]),(border,bmajticlocs[div])],\
    fill = (255,255,255,128))
    tick = str(bmajticks[div])
    d.text((border - len(tick)*5 - 40,bmajticlocs[div] - 10 ), \
    tick, font=fnt, fill=(255,255,255,128))
    
#axis labels
d.text((1.5*border + pixwt*0.5 - 150,pixhi + 2.5*border), \
"Time since ejection (Days)", font=fnt, fill=(255,255,255,128))
d.text((0.25*border - 10,0.75*border - 20), \
"Beta", font=fnt, fill=(255,255,255,128))

#plot title
plttitle = (comdenom.upper() + ' ' + comname[:-1] + ' dust tail')
tfnt = ImageFont.truetype(fontloc, 30)
d.text((1.5*border + pixwt*0.5 - len(plttitle)*5 - 100,.35*border), \
plttitle, font=tfnt, fill=(255,255,255,128))

dustimg.show()


