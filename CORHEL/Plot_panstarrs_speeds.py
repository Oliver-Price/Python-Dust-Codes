# -*- coding: utf-8 -*-
from astropy.time import Time
import easygui
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.ndimage

sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\FP Overplot")
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\Orbitplotting")

from orbitdata_loading_functions import *
from FP_plot_functions import *
from FP_diagnostics import *
from imagetime_methods import *
from conversion_routines import *
from io_methods import *
from simulation_setup import *
from particle_sim import *
from enlil_cooordinate_transform import xyz2rlatlon
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import astropy.units as u

t_plot = Time('2013-03-12T21:40',format='isot')
firstt = Time('2013-03-05T00:00',format='isot')
latt = Time('2013-03-21T11:50',format='isot')

#choosing comet data to use
inputfile = r'C:\\PhD\\Comet_data\\Input_files\\Input file_PANSTARRS_c2011l4_pt1.txt'

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

obsloc = 'Stereo_B'
 
#import the orbit data
comveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir)
eveceq = orb_vector_new(comdenom, 'Earth', pysav, orbitdir, opts = 'obs')

t_cell = abs(comveceq[:,0]-t_plot.jd).argmin()
e_cell = abs(comveceq[:,0]-latt.jd).argmin()

#%%

t_dust1 = Time('2013-03-05T12:00',format='isot')
beta1 = 0.8

t_dust2 = Time('2013-03-08T04:53',format='isot')
beta2 = 1.7

comcell1 = abs(comveceq[:,0] - t_dust1.jd).argmin()
comjds = Time(comveceq[:,0],format='jd').squeeze()

simt = 16*60*24 #20 days in minutes
nperday = 48
tstart1 = comveceq[comcell1,6:12]

(times1,positions1) = part_sim_fine_track(beta1, simt, nperday, tstart1)
jds1 = Time(t_dust1.jd+times1/1440,format='jd').squeeze()

comcell2 = abs(comveceq[:,0] - t_dust2.jd).argmin()

simt = 16*60*24 #20 days in minutes
nperday = 48
tstart2 = comveceq[comcell2,6:12]

(times2,positions2) = part_sim_fine_track(beta2, simt, nperday, tstart2)
jds2 = Time(t_dust2.jd+times2/1440,format='jd').squeeze()

#%%
r = np.linalg.norm(comveceq[:,6:9],axis=1)
times = Time(comveceq[:,0],format='jd')
timese = Time(comveceq[:,0],format='jd')

bpolc = np.genfromtxt ('bpol-2013-3-12-21-40-18.csv', delimiter=",")
vrc = np.genfromtxt ('vr-2013-3-12-21-40-18.csv', delimiter=",")

lonb = np.resize(bpolc[:,1],(59,180))
latb = np.resize(bpolc[:,2],(59,180))
polbr = np.resize(bpolc[:,3],(59,180))
polb = np.ones_like(polbr)
polb[polbr<0] = -1

lonv = np.resize(vrc[:,1],(59,180))
latv = np.resize(vrc[:,2],(59,180))
vr = np.resize(vrc[:,3],(59,180))

#%%

[f_r_hore,f_theta_hore,f_phi_hore] = xyz2rlatlon(comveceq[:,6],comveceq[:,7],comveceq[:,8])

f_hae = SkyCoord(f_phi_hore*u.deg, f_theta_hore*u.deg, f_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = times)
f_heeq = f_hae.transform_to(frames.HeliographicStonyhurst)

f_r_hnm = f_heeq.radius.value
f_phi_hnm = f_heeq.lon.value + 180
f_theta_hnm = f_heeq.lat.value

#%%

[e_r_hore,e_theta_hore,e_phi_hore] = xyz2rlatlon(eveceq[:,6],eveceq[:,7],eveceq[:,8])

e_hae = SkyCoord(e_phi_hore*u.deg, e_theta_hore*u.deg, e_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = timese)
e_heeq = e_hae.transform_to(frames.HeliographicStonyhurst)

e_r_hnm = e_heeq.radius.value
e_phi_hnm = e_heeq.lon.value + 180
e_theta_hnm = e_heeq.lat.value

#%%

[g_r_hore,g_theta_hore,g_phi_hore] = xyz2rlatlon(positions1[:,0],positions1[:,1],positions1[:,2])

g_hae = SkyCoord(g_phi_hore*u.deg, g_theta_hore*u.deg, g_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = jds1)
g_heeq = g_hae.transform_to(frames.HeliographicStonyhurst)

g_r_hnm = g_heeq.radius.value
g_phi_hnm = g_heeq.lon.value + 180
g_theta_hnm = g_heeq.lat.value

#%%

[p_r_hore,p_theta_hore,p_phi_hore] = xyz2rlatlon(positions2[:,0],positions2[:,1],positions2[:,2])

p_hae = SkyCoord(p_phi_hore*u.deg, p_theta_hore*u.deg, p_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = jds2)
p_heeq = p_hae.transform_to(frames.HeliographicStonyhurst)

p_r_hnm = p_heeq.radius.value
p_phi_hnm = p_heeq.lon.value + 180
p_theta_hnm = p_heeq.lat.value

#%%
fig, axs = plt.subplots(1, 1, figsize=(9,3))
test = axs.pcolormesh(lonv,latv,vr,shading='gouraud',cmap='plasma')
cbar = fig.colorbar(test, ax=axs)
axs.contour(lonb,latb,polb,colors='white',levels=0,linewidths=3)
goodlocs1 = abs(f_theta_hnm)<=60
goodlocs2 = times.jd>firstt.jd
goodlocs = goodlocs1&goodlocs2 
axs.plot(g_phi_hnm,g_theta_hnm,c='lime',lw=2)
axs.plot(p_phi_hnm,p_theta_hnm,c='magenta',lw=2)
axs.scatter(e_phi_hnm[goodlocs],e_theta_hnm[goodlocs],c='cyan',lw=2)
axs.plot(f_phi_hnm[goodlocs],f_theta_hnm[goodlocs],c='black',lw=3)
axs.set_ylim(-58,58)
axs.set_ylabel(r"Heliospheric Latitude (degrees)")
axs.set_xlabel(r"Heliospheric Longitude (degrees)")
cbar.set_label(r"Solar Wind Speed (km/s)", rotation=270,labelpad=15)
fig.subplots_adjust(left=0.09,right=1.05,bottom=0.2)

plt.savefig('C:\PhD\Python\Python-Dust-Codes\CORHEL\Panstarrs_VSW.png',dpi=300)