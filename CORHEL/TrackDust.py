# -*- coding: utf-8 -*-
from astropy.time import Time
import easygui
import os
import numpy as np
import sys
import matplotlib.pyplot as plt

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

t_dust = Time('2013-03-05T12:00',format='isot')
beta = 0.8

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
comveceq = orb_vector_new(comdenom, obsloc, pysav, orbitdir, opts = 'eq')

comcell = abs(comveceq[:,0] - t_dust.jd).argmin()

simt = 20*60*24 #10 days in minutes
nperday = 48
tstart = comveceq[comcell,6:12]

(times,positions) = part_sim_fine_track(beta, simt, nperday, tstart)
jds = Time(t_dust.jd+times/1440,format='jd').squeeze()

#%%
theta = np.radians(23.4392911)

ecliptic = np.zeros_like(positions[:,0:3])
ecliptic[:,0] = positions[:,0]
ecliptic[:,1] = positions[:,1]*np.cos(theta) + positions[:,2]*np.sin(theta)
ecliptic[:,2] = positions[:,2]*np.cos(theta) - positions[:,1]*np.sin(theta)

[f_r_hore,f_theta_hore,f_phi_hore] = xyz2rlatlon(ecliptic[:,0],ecliptic[:,1],ecliptic[:,2])

f_hae = SkyCoord(f_phi_hore*u.deg, f_theta_hore*u.deg, f_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = jds)
f_heeq = f_hae.transform_to(frames.HeliographicStonyhurst)

f_r_hnm = f_heeq.radius.value
f_phi_hnm = f_heeq.lon.value + 180
f_theta_hnm = f_heeq.lat.value

#%%
t_hcs =Time('2013-03-11T15:40',format='isot')
da = np.where((jds-t_hcs)>0)[0][0]
db = np.where((jds-t_hcs-0.25)>0)[0][0]

rr = 360/27.2753
dt = (jds[da:db]-t_hcs).jd
cor = rr*dt

t = np.array2string(jds[da:db].jd, formatter={'float_kind':lambda x: "%.8f" % x}, separator='\n', suppress_small=True)[1:-1]
a = np.array2string(f_r_hnm[da:db], formatter={'float_kind':lambda x: "%.8f" % x}, separator=',', suppress_small=True).replace("\n", "")[1:-1]
c = np.array2string(f_theta_hnm[da:db], formatter={'float_kind':lambda x: "%.8f" % x}, separator=',', suppress_small=True).replace("\n", "")[1:-1]
b = np.array2string(f_phi_hnm[da:db]-cor, formatter={'float_kind':lambda x: "%.8f" % x}, separator=',', suppress_small=True).replace("\n", "")[1:-1]

print(t_hcs.isot)
print("\n")
print(a)
print("\n")
print(b)
print("\n")
print(c)
print("\n")
print(t)