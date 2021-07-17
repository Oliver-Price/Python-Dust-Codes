# -*- coding: utf-8 -*-
import numpy as np
from astropy.time import Time
from enlil_cooordinate_transform import *
import numpy as np
import matplotlib.pyplot as plt

#%% FORWARD
csv = np.genfromtxt ('stereoa_hor.csv', delimiter=",")
[f_r_hore,f_theta_hore,f_phi_hore] = xyz2rlatlon(csv[:,1],csv[:,2],csv[:,3])

f_t = Time(csv[:,0],format='jd')
f_n = get_n(f_t)
f_L = get_L(f_n)
f_g = get_g(f_n)

f_phi_hae = f_phi_hore
f_theta_hae = f_theta_hore
f_r_hae = f_r_hore

f_lamb = get_lambda(f_L,f_g)
f_omega = get_omega(f_t)
f_bigtheta = get_theta(f_lamb,f_omega)

[f_x_hae,f_y_hae,f_z_hae] = rlatlon2xyz(f_r_hae,f_theta_hae,f_phi_hae)
[f_x_heeq,f_y_heeq,f_z_heeq] = hae2heeq(f_x_hae,f_y_hae,f_z_hae,f_bigtheta,f_omega)
[f_r_heeq,f_theta_heeq,f_phi_heeq] = xyz2rlatlon(f_x_heeq,f_y_heeq,f_z_heeq)
[f_r_hnm,f_theta_hnm,f_phi_hnm] = heeq2hnm(f_x_heeq,f_y_heeq,f_z_heeq)

f_sinbo = sind(f_lamb - f_omega) * sind(7.25)
f_cosbo = np.sqrt(1 - f_sinbo**2)
f_bo = np.degrees(np.arctan(f_sinbo/f_cosbo))

#%% BACK
csv2 = np.genfromtxt ('stereoa_rlatlon.csv', delimiter=",")
t_ref = Time('2013-02-22T03:03:00.099', format='isot', scale='utc') 
t_csv = t_ref.jd-0.00208+csv2[:,0]
b_t = Time(t_csv,format='jd')
b_n = get_n(b_t)
b_L = get_L(b_n)
b_g = get_g(b_n)
b_lamb = get_lambda(b_L,b_g)
b_omega = get_omega(b_t)
b_bigtheta = get_theta(b_lamb,b_omega)

b_r = csv2[:,1]
b_lat = csv2[:,2]
b_lon = csv2[:,3]

[b_x_heeq,b_y_heeq,b_z_heeq] = hnm2heeq(b_r,b_lat,b_lon)
[b_r_heeq,b_theta_heeq,b_phi_heeq] = xyz2rlatlon(b_x_heeq,b_y_heeq,b_z_heeq)
[b_x_hae,b_y_hae,b_z_hae] = heeq2hae(b_x_heeq,b_y_heeq,b_z_heeq,b_bigtheta,b_omega)
[b_r_hae,b_theta_hae,b_phi_hae] = xyz2rlatlon(b_x_hae,b_y_hae,b_z_hae)

b_sinbo = sind(b_lamb - b_omega) * sind(7.25)
b_cosbo = np.sqrt(1 - b_sinbo**2)
b_bo = np.degrees(np.arctan(b_sinbo/b_cosbo))

#%% verification

[b_x_hae1,b_y_hae1,b_z_hae1] =  rlatlon2xyz(b_r_hae,b_theta_hae,b_phi_hae)
[b_x_heeq1,b_y_heeq1,b_z_heeq1] = hae2heeq(b_x_hae1,b_y_hae1,b_z_hae1,b_bigtheta,b_omega)
[b_r1,b_lat1,b_lon1] = heeq2hnm(b_x_heeq1,b_y_heeq1,b_z_heeq1)

#%% plots

plt.plot(f_t.jd,f_r_hnm)
plt.plot(b_t.jd,b_r)