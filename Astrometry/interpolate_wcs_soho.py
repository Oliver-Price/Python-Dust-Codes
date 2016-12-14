import os
import numpy as np
from astropy.io import fits
from astropy import wcs
import datetime
import astropy.time

uncalsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Blue\processing'
solarxsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear\processing\Renamed Data'
radecsav = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Soho\C3_Clear'

uncal_list = sorted(os.listdir(uncalsav))
uncal_list = [s for s in uncal_list if ".fts" in s]
uncal_total = len(uncal_list)

cal_list = sorted(os.listdir(solarxsav))
cal_list = [s for s in cal_list if ".fits" in s]
cal_total = len(cal_list)

uncal_matches = [];uncal_tds = []

for i in xrange(0,cal_total):
    csec = int(cal_list[i][13:15])
    cmin = int(cal_list[i][11:13])
    chour = int(cal_list[i][9:11])
    cday = int(cal_list[i][6:8])
    cmonth = int(cal_list[i][4:6])
    cyear = int(cal_list[i][0:4])
    ctime = astropy.time.Time(datetime.datetime(cyear,cmonth,cday,chour,cmin,csec))
    td = 1000
    for j in xrange(0,uncal_total):
        usec = int(uncal_list[j][11:13])
        umin = int(uncal_list[j][9:11])
        uhour = int(uncal_list[j][7:9])
        uday = int(uncal_list[j][4:6])
        umonth = int(uncal_list[j][2:4])
        uyear = 2000 + int(uncal_list[j][0:2])
        utime = astropy.time.Time(datetime.datetime(uyear,umonth,uday,uhour,umin,usec))
        if abs(ctime.jd - utime.jd) < td:
            td = abs(ctime.jd - utime.jd)
            match = uncal_list[j]
    uncal_matches.append(match)
    uncal_tds.append(td*24*60)

#%%

sx_fits = os.path.join(solarxsav,cal_list[0])
rd_fits = os.path.join(radecsav,cal_list[0])
uc_fits = os.path.join(uncalsav,uncal_list[0])

sx_img = fits.open(sx_fits)
rd_img = fits.open(rd_fits)
uc_img = fits.open(uc_fits)

ws = wcs.WCS(sx_img[0].header)
wr = wcs.WCS(rd_img[0].header)
wu = wcs.WCS(uc_img[0].header)

#make a 2xN array of all pixel locations
ya = rd_img[0].data.shape[0]
xa = rd_img[0].data.shape[1]
xv, yv = np.meshgrid(np.arange(xa), np.arange(ya))
xv = np.reshape(xv, (xa*ya,1))
yv = np.reshape(yv, (xa*ya,1))
coords = np.zeros((xa*ya,2))
coords[:,0] = xv[:,0]
coords[:,1] = yv[:,0]

#convert each pixel locaion to an RA and DEC array
radecs = wr.wcs_pix2world(coords, 0)
ra = np.reshape(radecs[:,0], (ya,xa))
dec = np.reshape(radecs[:,1], (ya,xa))

#make a 2xN array of all pixel locations
uya = uc_img[0].data.shape[0]
uxa = uc_img[0].data.shape[1]
uxv, uyv = np.meshgrid(np.arange(uxa), np.arange(uya))
uxv = np.reshape(uxv, (uxa*uya,1))
uyv = np.reshape(uyv, (uxa*uya,1))
ucoords = np.zeros((uxa*uya,2))
ucoords[:,0] = uxv[:,0]
ucoords[:,1] = uyv[:,0]

solxys = ws.wcs_pix2world(coords, 0)
solar_x = np.reshape(solxys[:,0], (ya,xa))
solar_x[np.where(solar_x > 180)] = solar_x[np.where(solar_x > 180)] - 360
solar_y = np.reshape(solxys[:,1], (ya,xa))

solar_xa = solar_x*3600
solar_ya = solar_y*3600

#convert each pixel locaion to an RA and DEC array
uc_xy = wu.wcs_pix2world(ucoords, 0)
u_solar_rawx = -np.reshape(uc_xy[:,0], (uya,uxa))
u_solar_rawy = -np.reshape(uc_xy[:,1], (uya,uxa))

u_solar_x = (360 + u_solar_rawx/3600)%360
u_solar_y = u_solar_rawy/3600
u_solay_xy = np.zeros((uxa*uya,2), dtype=float)
u_solay_xy[:,0] = np.reshape(u_solar_x, (uxa*uya,1))[:,0]
u_solay_xy[:,1]  = np.reshape(u_solar_y, (uxa*uya,1))[:,0]

pixlocs = ws.wcs_world2pix(u_solay_xy,0)
pixlocs[:,1] = abs(pixlocs[:,1] - 1024) + 1

u_radecs = wr.wcs_pix2world(pixlocs, 0)
u_ra = np.reshape(u_radecs[:,0], (uya,uxa))
u_dec = np.reshape(u_radecs[:,1], (uya,uxa))

pix_x = np.reshape(pixlocs[:,0], (uya,uxa))
pix_y = np.reshape(pixlocs[:,1], (uya,uxa))
