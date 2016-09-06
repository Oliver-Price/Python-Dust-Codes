#auxiliary conversion functions
import math as m
import numpy as np
import easygui
from astropy.io import fits
import os
import sys

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

#%% - fix occurences where range is incorrect or where data wraps from 360 ra to 0

def fixwraps(ra, ramax, ramin): 
    ra_m = np.copy(ra%360)
    if (ramax-ramin) > 270:
        circlocs = np.where(ra_m < 180)
        ra_m[circlocs] = ra_m[circlocs] + 360
        rafmin = np.amin(ra_m)
        rafmax = np.amax(ra_m)
    else:
        ra_m = ra
        rafmin = ramin
        rafmax = ramax
    
    return ra_m, rafmin, rafmax
 
#%% - ensure data input is in correct form
   
def correct_for_imagetype(imagedir, fitsin, fitsinfile):
    #choose img type
    colormsg = "Please select image type:"
    colorchoices = ["Colour","Black & White"]
    reply = easygui.buttonbox(colormsg, choices=colorchoices)
        
    if reply == "Colour":

        #we have to use a temporary 1d fits file to get WCS data
        #this is due to astropy being unable to handle RGB fits images well
        fitstemp = os.path.join(imagedir, 'temporary_' + fitsinfile)
        hdulist = fits.open(fitsin)
        
        #opening color data
        colours = (hdulist[0].data)
        colr = colours[0,:,:]
        colg = colours[1,:,:]
        colb = colours[2,:,:]
        
        #making a 1d image for doing ra/dec coordinates
        oldtemp = os.path.isfile(fitstemp) #if one doesn't already exist
        if oldtemp == False: 
            hdulist[0].data = hdulist[0].data[0] #makes image from red plane
            hdulist[0].header['naxis'] = 1       #edits header so ds9 will work
            hdulist.writeto(fitstemp)
            
        hdulist.close()
        fitscoords = fitstemp #directs program to look at temp image for coords
     
    elif reply == "Black & White":
        
        #simple case, we can use base image for coords
        hdulist = fits.open(fitsin)
        colours = (hdulist[0].data)
        fitscoords = fitsin
        colr = colours #black and white colour scheme
        colg = colours
        colb = colours
        
        return colr, colg, colb, fitscoords
 
#%% corrects for data > 255, inversions, wrong data type etc.
       
def col_corrections(inv,colr,colg,colb):
    
    im_max = max(np.max(colr),np.max(colg),np.max(colb))
    if im_max > 255:   
        coltr = np.rint(colr*255/im_max)
        coltg = np.rint(colg*255/im_max)
        coltb = np.rint(colb*255/im_max)
    else:
        coltr = colr; coltg = colg; coltb = colb
        
    if inv == False:   
        colcr = coltr; colcg = coltg; colcb = coltb
        backgr_fill = (0,0,0,255)
        featur_fill = (255,255,255,255)
    elif inv == True:
        colcr = 255 - coltr; colcg = 255  - coltg; colcb = 255 - coltb
        backgr_fill = (255,255,255,255)
        featur_fill = (0,0,0,255)
        
    if 'float' in str(colr.dtype):
         colcr = colcr.astype(int)
         colcg = colcg.astype(int)
         colcb = colcb.astype(int)
         
    return backgr_fill, featur_fill, colcr, colcg, colcb

#%% load in correct observer location
    
def get_obs_loc(obslocstr, imagedir):
    
    bool_locs = np.array([(('EARTH' in obslocstr) or ('Earth' in obslocstr)),
                 (('STEREO-A' in obslocstr) or ('Stereo-A' in obslocstr)
                 or ('STEREO_A' in obslocstr) or ('Stereo_A' in obslocstr)),
                 (('STEREO-B' in obslocstr) or ('Stereo-B' in obslocstr)
                 or ('STEREO_B' in obslocstr) or ('Stereo_B' in obslocstr)),
                 (('SOHO' in obslocstr) or ('Soho' in obslocstr))])
    name_locs = np.array(['Earth', 'Stereo_A', 'Stereo_B', 'Soho'])
    case_locs = np.size(np.nonzero(bool_locs))
    if case_locs > 1:
        obsmsg = "Please select observer location"
        obschoices = name_locs[bool_locs].tolist()
        obsloc = easygui.buttonbox(obsmsg, choices=obschoices)
        imagedir = os.path.join(imagedir, obsloc)
    elif case_locs == 1:
        obsloc = name_locs[bool_locs][0]
    else: sys.exit("No Good Observer Location")
    
    return obsloc, imagedir
    
#%% load in correct observer location
    
def get_stereo_instrument(imagedir):

    name_locs = np.array(['HI-1', 'HI-2'])
    stermsg = "Please select STEREO instrument"
    sterchoices = name_locs.tolist()
    sterinst = easygui.buttonbox(stermsg, choices=sterchoices)
    imagedir = os.path.join(imagedir, sterinst)

    return sterinst, imagedir
    
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
        return start_locs[0][l_index], stop_locs[0][l_index]
    else: return None, None
        