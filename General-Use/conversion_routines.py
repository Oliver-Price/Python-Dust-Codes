#auxiliary conversion functions
import math as m
import numpy as np
import easygui
import os
import sys
from math import log10, floor

#%% - converts month text to a number

def mon2num(month):
    if month == 'Jan':
        mnum = 1
    elif month == 'Feb':
        mnum = 2
    elif month == 'Mar':
        mnum = 3
    elif month == 'Apr':
        mnum = 4
    elif month == 'May':
        mnum = 5
    elif month == 'Jun':
        mnum = 6
    elif month == 'Jul':
        mnum = 7
    elif month == 'Aug':
        mnum = 8
    elif month == 'Sep':
        mnum = 9
    elif month == 'Oct':
        mnum = 10
    elif month == 'Nov':
        mnum = 11
    elif month == 'Dec':
        mnum = 12
    return mnum
    
#%% - converts earth centric equatorial xyz coord to ra and dec
    
def pos2radec(position,fourpi = False):
    ra = m.atan(position[1]/position[0])*360/(2*m.pi)
    if position[0] >= 0:
        if position[1] < 0:
            ra = ra + 360
    if position[0] < 0:
        ra = ra + 180
    r = np.linalg.norm(position)
    dec = m.asin(position[2]/r)*360/(2*m.pi)
    if fourpi == False: return (ra,dec)
    elif fourpi == True: 
        if ra < 180: return (ra + 360,dec)
        if ra >= 180: return (ra,dec)

#%% - fix occurences where range is incorrect or where data wraps from 360 ra to 0

def fixwraps(ra, ramax, ramin): 
    ra_m = np.copy(ra%360)
    if (ramax-ramin) > 270:
        circlocs = np.where(ra_m < 180)
        ra_m[circlocs] = ra_m[circlocs] + 360
    else: pass
    rafmin = np.amin(ra_m)
    rafmax = np.amax(ra_m)
    return ra_m, rafmin, rafmax
    
#%% exactly what it says
   
def find_largest_nonzero_block(array_in):
    
    fake_edge = np.array([0])
    array_in = np.concatenate((fake_edge, array_in))
    array_in = np.concatenate((array_in, fake_edge))
    
    diff_arr = np.diff(array_in)
    start_locs = np.where(diff_arr == 1)
    stop_locs = np.where(diff_arr == -1)
    
    no_blocs = np.size(stop_locs[0])
    sizes = np.empty((no_blocs),dtype = int)
    for x in xrange(0,no_blocs):
        sizes[x] = stop_locs[0][x] - start_locs[0][x]
    
    if sizes.size > 0:
        l_index = np.argmax(sizes)
        return start_locs[0][l_index], stop_locs[0][l_index]-1
    else: return None, None

#%% rounding functions wot i copied from stack overflow

#to 1 significant figure
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

#to nearest multiple of base
def round_to_base(x, base):
    return base * round(float(x)/base)
