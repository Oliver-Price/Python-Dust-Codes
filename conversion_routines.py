#auxiliary conversion functions
import math as m
import numpy as np

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
    
def pos2radec(position):
    ra = m.atan(position[1]/position[0])*360/(2*m.pi)
    if position[0] >= 0:
        if position[1] < 0:
            ra = ra + 360
    if position[0] < 0:
        ra = ra + 180
    r = np.linalg.norm(position)
    dec = m.asin(position[2]/r)*360/(2*m.pi)
    return (ra,dec)
    
    
    