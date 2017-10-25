# -*- coding: utf-8 -*-
from astropy.io import fits
from astropy import wcs
import astropy.time
import os
import numpy as np
import sys
import datetime
import copy

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")

from orbitdata_loading_functions import orb_vector, orb_obs
from imagetime_methods import image_time_stereo
from conversion_routines import pos2radec

#%%

#uncalsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue\processing'
uncalsav = r"C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue\processing\Updated\Renamed"
calsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear'
calbluefinal = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue\processing\Updated\Calibrated'

uncal_list = sorted(os.listdir(uncalsav))
uncal_list = [s for s in uncal_list if ".f" in s]
uncal_total = len(uncal_list)

cal_list = sorted(os.listdir(calsav))
cal_list = [s for s in cal_list if ".f" in s]
cal_total = len(cal_list)

uncal_matches = [];tds = []; cal_matches = []

for j in range(0,uncal_total):
    usec = int(uncal_list[j][13:15])
    umin = int(uncal_list[j][11:13])
    uhour = int(uncal_list[j][9:11])
    uday = int(uncal_list[j][6:8])
    umonth = int(uncal_list[j][4:6])
    uyear = int(uncal_list[j][0:4])
    utime = astropy.time.Time(datetime.datetime(uyear,umonth,uday,uhour,umin,usec))
    td = 1000
    for i in range(0,cal_total):
        csec = int(cal_list[i][13:15])
        cmin = int(cal_list[i][11:13])
        chour = int(cal_list[i][9:11])
        cday = int(cal_list[i][6:8])
        cmonth = int(cal_list[i][4:6])
        cyear = int(cal_list[i][0:4])
        ctime = astropy.time.Time(datetime.datetime(cyear,cmonth,cday,chour,cmin,csec))
        if abs(ctime.jd - utime.jd) < td:
            td = abs(ctime.jd - utime.jd)
            calmatch = cal_list[i]
            uncalmatch = uncal_list[j]
    if td*24*60 < 15:
        uncal_matches.append(uncalmatch)
        cal_matches.append(calmatch)
        tds.append(td*24*60)

#%%
#choosing comet data to use
inputfile = "C:\PhD\Comet_data\Input_files\Input file_McNaught_c2006p1_pt1.txt"
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

#choose observer locations
obsloc = 'Soho'

#import the orbit data
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,eq')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'eq')
comveceq10 = orb_vector(comdenom, obsloc, pysav, orbitdir,
                        horiztag, opts = 'eq,d10')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag)

#%%

for x in range(len(cal_matches)):

    [ctime,uncertainty_range_exists] = image_time_stereo(cal_matches[x])
    
    comcel = np.where(abs(obsveceq[:,0] - ctime.jd)==abs(obsveceq[:,0] - ctime.jd).min())[0][0]
    sun_pos = pos2radec(-obsveceq[comcel,6:9],False)
    
    cal_soho = os.path.join(calsav,cal_matches[x])
    uncal_soho = os.path.join(uncalsav,uncal_matches[x])
    
    data_cal, header_cal = fits.getdata(cal_soho, header=True)
    data_un, header_un = fits.getdata(uncal_soho , header=True)
    
    w = wcs.WCS(header_cal)
    sun_coords = w.wcs_world2pix([sun_pos],0)
    
    header_nu = copy.copy(header_cal)
    
    header_nu['BITPIX'] = header_un['BITPIX']
    header_nu['NAXIS1'] = header_un['NAXIS1']
    header_nu['NAXIS2'] = header_un['NAXIS2']
    header_nu['DATE-OBS'] = header_un['DATE-OBS']
    header_nu['TIME-OBS'] = header_un['TIME-OBS']
    header_nu['CRPIX1'] = header_un['NAXIS1'] - header_un['CRPIX1'] + header_cal['CRPIX1'] - sun_coords[0][0]
    header_nu['CRPIX2'] = header_un['NAXIS2'] - header_un['CRPIX2'] + header_cal['CRPIX2'] - sun_coords[0][1]
    
    blue_sav = os.path.join(calbluefinal, uncal_matches[x])
    fits.writeto(blue_sav, np.flipud(np.fliplr(data_un)), header_nu, clobber=True)
    
#%%
'''
exptimes = np.zeros((len(uncal_matches)),dtype=float)
for fits_no in range(0,len(uncal_matches)):
    
    fits_wcs_loc = os.path.join(uncalsav, uncal_matches[fits_no])
    
    data_wcs, header_wcs = fits.getdata(fits_wcs_loc, header=True)

    exptimes[fits_no] = header_wcs['EXPTIME']
'''