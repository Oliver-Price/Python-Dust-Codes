# -*- coding: utf-8 -*-
'''
Coordinate transforms for ENLIL

Earth at:
'''
import numpy as np
import math as m
import astropy.time
import datetime

def cosd(p):
    return np.cos(np.radians(p))

def sind(q):
    return np.sin(np.radians(q))

#%% HEEQ

#theta = lat, phi = long
def hnm2heeq(r,theta,phi):
    x = r * sind(theta) * cosd(phi + 180)
    y = r * sind(theta) * sind(phi + 180)
    z = r * cosd(theta)
    return (x,y,z)

#theta = lat, phi = long
def heeq2heeqxyz(r,theta,phi):
    thetau = 90 - theta
    x = r * sind(thetau) * cosd(phi)
    y = r * sind(thetau) * sind(phi)
    z = r * cosd(thetau)
    return (x,y,z)

#%% SUN ORBITAL PARS

cos_i = cosd(7.25)
sin_i = sind(7.25)

def get_n(time):
    return (time.jd - 2451545.0)

def get_L(n):
    return (280.460 + 0.9856474*n)%360

def get_g(n):
    return (357.528 + 0.9856003*n)%360

def get_lambda(L,g):
    return (L + 1.915*sind(g) + 0.02*sind(2*g))%360

def get_omega(time):
    return (73.6667 + 0.013958*(time.mjd + 3243.75)/365.25)%360

def get_theta(lb,o):
    return np.degrees( m.atan(cos_i * m.tan( m.radians(lb - o) ) ) )%360

#%% HAE

def heeq2hae(a,b,c,theta,omega):
    x = a*(cosd(theta)*cosd(omega) - sind(theta)*sind(omega)*cos_i)
    x += b*(-cosd(theta)*sind(omega) - sind(theta)*cosd(omega)*cos_i)
    x += c*sin_i*sind(omega)
    y = a*(sind(theta)*cosd(omega) + cosd(theta)*sind(omega)*cos_i)
    y += b*(-sind(theta)*sind(omega) + cosd(theta)*cosd(omega)*cos_i)
    y += -c*sin_i*cosd(omega)
    z = a*sin_i*sind(theta)
    z += b*sin_i*cosd(theta)
    z += c*cos_i
    return (x,y,z)