# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 17:54:46 2018

@author: Ollie
"""
import numpy as np

i = -23.43929
z = 115.912
x = 81.706

earthframe = comveceq

cosi = np.cos(np.radians(i))
sini = np.sin(np.radians(i))
cosz = np.cos(np.radians(z))
sinz = np.sin(np.radians(z))
cosx = np.cos(np.radians(x))
sinx = np.sin(np.radians(x))

cosi2 = np.cos(np.radians(-i))
sini2 = np.sin(np.radians(-i))
cosz2 = np.cos(np.radians(-z))
sinz2 = np.sin(np.radians(-z))
cosx2 = np.cos(np.radians(-x))
sinx2 = np.sin(np.radians(-x))

temp = np.empty_like(earthframe[...,6:12])
temp2 = np.empty_like(earthframe[...,6:12])
temp5 = np.empty_like(earthframe[...,6:12])
temp6 = np.empty_like(earthframe[...,6:12])
orb_plane = np.empty_like(earthframe[...,6:12])
eframe = np.empty_like(earthframe[...,6:12])

#rotate about i (from earth's frame to ecliptic) = rotation by earth's tilt
temp[...,0] = earthframe[...,6]
temp[...,1] = earthframe[...,7]*cosi - earthframe[...,8]*sini
temp[...,2] = earthframe[...,7]*sini + earthframe[...,8]*cosi
temp[...,3] = earthframe[...,9]
temp[...,4] = earthframe[...,10]*cosi - earthframe[...,11]*sini
temp[...,5] = earthframe[...,10]*sini + earthframe[...,11]*cosi

#rotate about z (up from ecliptic) = rotation by ascending node (may be off given value by 180)
temp2[...,0] = temp[...,0]*cosz - temp[...,1]*sinz
temp2[...,1] = temp[...,0]*sinz + temp[...,1]*cosz
temp2[...,2] = temp[...,2]
temp2[...,3] = temp[...,3]*cosz - temp[...,4]*sinz
temp2[...,4] = temp[...,3]*sinz + temp[...,4]*cosz
temp2[...,5] = temp[...,5]

#rotate about x (pointing at ascending node) = rotation by inclination
orb_plane[...,0] = temp2[...,0]
orb_plane[...,1] = temp2[...,1]*cosx - temp2[...,2]*sinx
orb_plane[...,2] = temp2[...,1]*sinx + temp2[...,2]*cosx
orb_plane[...,3] = temp2[...,3]
orb_plane[...,4] = temp2[...,4]*cosx - temp2[...,5]*sinx
orb_plane[...,5] = temp2[...,4]*sinx + temp2[...,5]*cosx

orbframe = orb_plane

#rotate about x (pointing at ascending node) = rotation by inclination
temp6[...,0] = orbframe[...,0]
temp6[...,1] = orbframe[...,1]*cosx2 - orbframe[...,2]*sinx2
temp6[...,2] = orbframe[...,1]*sinx2 + orbframe[...,2]*cosx2
temp6[...,3] = orbframe[...,3]
temp6[...,4] = orbframe[...,4]*cosx2 - orbframe[...,5]*sinx2
temp6[...,5] = orbframe[...,4]*sinx2 + orbframe[...,5]*cosx2

#rotate about z (up from ecliptic) = rotation by ascending node (may be off given value by 180)
temp5[...,0] = temp6[...,0]*cosz2 - temp6[...,1]*sinz2
temp5[...,1] = temp6[...,0]*sinz2 + temp6[...,1]*cosz2
temp5[...,2] = temp6[...,2]
temp5[...,3] = temp6[...,3]*cosz2 - temp6[...,4]*sinz2
temp5[...,4] = temp6[...,3]*sinz2 + temp6[...,4]*cosz2
temp5[...,5] = temp6[...,5]

#rotate about i (from earth's frame to ecliptic) = rotation by earth's tilt
eframe[...,0] = temp5[...,0]
eframe[...,1] = temp5[...,1]*cosi2 - temp5[...,2]*sini2
eframe[...,2] = temp5[...,1]*sini2 + temp5[...,2]*cosi2
eframe[...,3] = temp5[...,3]
eframe[...,4] = temp5[...,4]*cosi2 - temp5[...,5]*sini2
eframe[...,5] = temp5[...,4]*sini2 + temp5[...,5]*cosi2

out = earthframe[...,6:12] - eframe