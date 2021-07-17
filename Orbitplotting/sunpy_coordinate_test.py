# -*- coding: utf-8 -*-
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
from enlil_cooordinate_transform import xyz2rlatlon
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

#%% FORWARD
csv = np.genfromtxt ('earth_hor.csv', delimiter=",")
[f_r_hore,f_theta_hore,f_phi_hore] = xyz2rlatlon(csv[:,1],csv[:,2],csv[:,3])
f_t = Time(csv[:,0],format='jd')

f_hae = SkyCoord(f_phi_hore*u.deg, f_theta_hore*u.deg, f_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = f_t)
f_heeq = f_hae.transform_to(frames.HeliographicStonyhurst)

#%% BACK
csv2 = np.genfromtxt ('earth_rlatlon.csv', delimiter=",")
t_ref = Time('2013-02-22T03:03:00.099', format='isot', scale='utc') 
t_csv = t_ref.jd-0.00208+csv2[:,0]
b_t = Time(t_csv,format='jd')

b_r = csv2[:,1]
b_lat = csv2[:,2]
b_lon = csv2[:,3] - 180

#%%
plt.plot(f_t.jd,f_heeq.lat.value)
plt.plot(b_t.jd,b_lat)

#%%
plt.plot(f_t.jd,f_heeq.lon.value)
plt.plot(b_t.jd,b_lon)

#%%
plt.plot(f_t.jd,f_heeq.radius.value)
plt.plot(b_t.jd,b_r)

#%% FORWARD
csv = np.genfromtxt ('stereoa_hor.csv', delimiter=",")
[f_r_hore,f_theta_hore,f_phi_hore] = xyz2rlatlon(csv[:,1],csv[:,2],csv[:,3])
f_t = Time(csv[:,0],format='jd')

f_hae = SkyCoord(f_phi_hore*u.deg, f_theta_hore*u.deg, f_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = f_t)
f_heeq = f_hae.transform_to(frames.HeliographicStonyhurst)

#%% BACK
csv2 = np.genfromtxt ('stereoa_rlatlon.csv', delimiter=",")
t_ref = Time('2013-02-22T03:03:00.099', format='isot', scale='utc') 
t_csv = t_ref.jd-0.00208+csv2[:,0]
b_t = Time(t_csv,format='jd')

b_r = csv2[:,1]
b_lat = csv2[:,2]
b_lon = csv2[:,3] - 180

#%%
plt.plot(f_t.jd,f_heeq.lat.value)
plt.plot(b_t.jd,b_lat)

#%%
plt.plot(f_t.jd,f_heeq.lon.value)
plt.plot(b_t.jd,b_lon)

#%%
plt.plot(f_t.jd,f_heeq.radius.value)
plt.plot(b_t.jd,b_r)