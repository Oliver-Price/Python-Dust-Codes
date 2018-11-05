# -*- coding: utf-8 -*-

import numpy as np
def earthplane2orbplane(x,z,earthframe):

	i = -23.43929
	
	cosi = np.cos(np.radians(i))
	sini = np.sin(np.radians(i))
	cosz = np.cos(np.radians(z))
	sinz = np.sin(np.radians(z))
	cosx = np.cos(np.radians(x))
	sinx = np.sin(np.radians(x))

	temp = np.empty_like(earthframe[...,6:12])
	temp2 = np.empty_like(earthframe[...,6:12])
	orb_plane = np.empty_like(earthframe[...,6:12])

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
	
	return orb_plane
	
def orbplane2earthplane(x,z,orbframe):

	i = 23.43929

	cosi = np.cos(np.radians(i))
	sini = np.sin(np.radians(i))
	cosz = np.cos(np.radians(z))
	sinz = np.sin(np.radians(z))
	cosx = np.cos(np.radians(x))
	sinx = np.sin(np.radians(x))

	temp = np.empty_like(orbframe)
	temp2 = np.empty_like(orbframe)
	eframe = np.empty_like(orbframe)
	
	#rotate about x (pointing at ascending node) = rotation by inclination
	temp[...,0] = orbframe[...,0]
	temp[...,1] = orbframe[...,1]*cosx - orbframe[...,2]*sinx
	temp[...,2] = orbframe[...,1]*sinx + orbframe[...,2]*cosx
	temp[...,3] = orbframe[...,3]
	temp[...,4] = orbframe[...,4]*cosx - orbframe[...,5]*sinx
	temp[...,5] = orbframe[...,4]*sinx + orbframe[...,5]*cosx

	#rotate about z (up from ecliptic) = rotation by ascending node (may be off given value by 180)
	temp2[...,0] = temp[...,0]*cosz - temp[...,1]*sinz
	temp2[...,1] = temp[...,0]*sinz + temp[...,1]*cosz
	temp2[...,2] = temp[...,2]
	temp2[...,3] = temp[...,3]*cosz - temp[...,4]*sinz
	temp2[...,4] = temp[...,3]*sinz + temp[...,4]*cosz
	temp2[...,5] = temp[...,5]

	#rotate about i (from earth's frame to ecliptic) = rotation by earth's tilt
	eframe[...,0] = temp2[...,0]
	eframe[...,1] = temp2[...,1]*cosi - temp2[...,2]*sini
	eframe[...,2] = temp2[...,1]*sini + temp2[...,2]*cosi
	eframe[...,3] = temp2[...,3]
	eframe[...,4] = temp2[...,4]*cosi - temp2[...,5]*sini
	eframe[...,5] = temp2[...,4]*sini + temp2[...,5]*cosi
	
	return eframe
