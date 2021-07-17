# -*- coding: utf-8 -*-
import os
os.chdir("C:\PhD\Python\Python-Dust-Codes\Orbitplotting")
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
from enlil_cooordinate_transform import xyz2rlatlon
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames


csv = np.genfromtxt ('panstarrs_hor.csv', delimiter=",")
[f_r_hore,f_theta_hore,f_phi_hore] = xyz2rlatlon(csv[:,1],csv[:,2],csv[:,3])
f_t = Time(csv[:,0],format='jd')

f_hae = SkyCoord(f_phi_hore*u.deg, f_theta_hore*u.deg, f_r_hore*u.au, frame='heliocentricmeanecliptic', obstime = f_t)
f_heeq = f_hae.transform_to(frames.HeliographicStonyhurst)

f_r_hnm = f_heeq.radius.value
f_phi_hnm = f_heeq.lon.value + 180
f_theta_hnm = f_heeq.lat.value

#%%
plt.plot(f_t.jd[408:552],f_r_hnm[408:552])

#%%
plt.plot(f_t.jd[408:552],f_phi_hnm[408:552])

#%%
plt.plot(f_t.jd[408:552],f_theta_hnm[408:552])

#%%
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(csv[:,1],csv[:,2],csv[:,3])

#%%

t_hcs =Time('2013-03-10T21:40',format='isot')
da = np.argmin(np.abs(f_t-t_hcs))
dh = 6

rr = 360/27.2753
dt = (f_t[da:da+6]-t_hcs).jd
cor = rr*dt

t = np.array2string(f_t[da:da+6].jd, formatter={'float_kind':lambda x: "%.8f" % x}, separator='\n', suppress_small=True)[1:-1]
a = np.array2string(f_r_hnm[da:da+dh], formatter={'float_kind':lambda x: "%.8f" % x}, separator=',', suppress_small=True).replace("\n", "")[1:-1]
c = np.array2string(f_theta_hnm[da:da+dh], formatter={'float_kind':lambda x: "%.8f" % x}, separator=',', suppress_small=True).replace("\n", "")[1:-1]
b = np.array2string(f_phi_hnm[da:da+dh]-cor, formatter={'float_kind':lambda x: "%.8f" % x}, separator=',', suppress_small=True).replace("\n", "")[1:-1]

print(t_hcs.isot)
print(a)
print(b)
print(c)
print(t)