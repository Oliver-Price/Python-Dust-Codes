# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from astropy.io import fits
import datetime
import astropy.time

betau = 2.5
betal = 0.1
bno = 400
simtu = 6.0
simtl = 1.0
tno = 400

datafolder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Linux Data\SOHO Blue\fits'
plotsave = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Linux Data\SOHO Blue\plots'
image_list = sorted(os.listdir(datafolder))
image_total = len(image_list)

for image_id in range(0,image_total):

    image_procname = image_list[image_id].split('.')[0][:-7]
    image_mapped = os.path.join(datafolder,image_list[image_id])
    print(image_procname)
    
    hdulist = fits.open(image_mapped)
    srcolors = (hdulist[0].data).astype(float)
    srcolors = np.swapaxes(srcolors,0,1)
    
    tvals = np.linspace(simtl, simtu, tno)
    bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))
    
    csec = int(image_procname[13:15])
    cmin = int(image_procname[11:13])
    chour = int(image_procname[9:11])
    cday = int(image_procname[6:8])
    cmonth = int(image_procname[4:6])
    cyear = int(image_procname[0:4])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                chour , cmin, csec))
    
    #%%
    
    betaval = 0.7
    
    row = abs(bvals-betaval).argmin()
    rowS = np.where(srcolors[:,row]<=5000)[0][0]
    rowE = np.where(srcolors[:,row]>=400)[0][-1]
    
    betaval2 = 1.6
    
    row2 = abs(bvals-betaval2).argmin()
    rowS2 = np.where(srcolors[:,row2]<=5000)[0][0]
    rowE2 = np.where(srcolors[:,row2]>=400)[0][-1]
    
    #%%
    
    fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1)
    plottitle = image_procname[:15]
    
    plt.tight_layout(h_pad=5.0,rect=[0, 0.03, 1, 0.95])
    plt.subplots_adjust(top=0.85)
    plt.suptitle(plottitle, fontsize=16)
    ax1.semilogy(ctime.jd - tvals[rowS:rowE],srcolors[rowS:rowE,row])
    ax1.set_title("beta = " + str(betaval))
    ax1.set_xlabel('Julian Date')
    ax1.set_ylabel('Brightness')
    ax1.set_ylim(400,5000)
    ax1.set_xlim(0.6+2.45411e6,5.0+2.45411e6)
    ax1.axvline(x=(ctime.jd - 2.5),c='g')
    
    ax2.semilogy(ctime.jd - tvals[rowS2:rowE2],srcolors[rowS2:rowE2,row2])
    ax2.set_title("beta = " + str(betaval2))
    ax2.set_xlabel('Julian Date')
    ax2.set_ylabel('Brightness')
    ax2.set_ylim(400,5000)
    ax2.set_xlim(1.2+2.45411e6,5.2+2.45411e6)
    ax2.axvline(x=(ctime.jd - 2.1),c='r')
        
    #plt.show()
    fig.subplots_adjust(bottom=0.15,left=0.15)
    
    plotsavename = os.path.join(plotsave,image_procname[:15] + '.png')
    fig.savefig(plotsavename,dpi=150)
