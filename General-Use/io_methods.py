#auxiliary functions for easy IO
import easygui
from astropy.io import fits
import os
import numpy as np
import sys
 
#%% - methods deals with various size data cubes
  
def correct_for_imagetype(imagedir, fitsin, fitsinfile):

    hdulist = fits.open(fitsin)
    colours = (hdulist[0].data)
    no_dim = colours.ndim
    
    if no_dim == 2:
        fitscoords = fitsin
        colr = colours #black and white colour scheme
        colg = colours
        colb = colours
    
    elif no_dim == 3:
        
        colr = colours[0,:,:]
        colg = colours[1,:,:]
        colb = colours[2,:,:]
    
        #making a 1d image for doing ra/dec coordinates
        fitstemp = os.path.join(imagedir, 'temporary_' + fitsinfile)
        oldtemp = os.path.isfile(fitstemp) #if one doesn't already exist
        if oldtemp == False: 
            hdulist[0].data = hdulist[0].data[0] #makes image from red plane
            hdulist[0].header['naxis'] = 1       #edits header so ds9 will work
            hdulist.writeto(fitstemp)

        fitscoords = fitstemp #directs program to look at temp image for coords
    
    return colr, colg, colb, fitscoords
        
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

    sname_locs = np.array(['HI-1', 'HI-2'])
    stermsg = "Please select STEREO instrument"
    sterchoices = sname_locs.tolist()
    sterinst = easygui.buttonbox(stermsg, choices=sterchoices)
    
    mname_locs = np.array(['Base', 'Difference', 'Multiscale Gaussian Normalised'])
    modemsg = "Please select image type"
    modechoices = mname_locs.tolist()
    mode = easygui.buttonbox(modemsg, choices=modechoices)
    
    if "Diff" in mode: sterinst = sterinst + '-diff'
    if "Multi" in mode: sterinst = sterinst + '-MGN'
    imagedir = os.path.join(imagedir, sterinst)

    return sterinst, imagedir
    
def get_soho_instrument(imagedir):

    sohoinst = 'C3_Clear'
    
    mname_locs = np.array(['Base', 'Difference', 'Multiscale Gaussian Normalised'])
    modemsg = "Please select image type"
    modechoices = mname_locs.tolist()
    mode = easygui.buttonbox(modemsg, choices=modechoices)
    
    if "Diff" in mode: sohoinst = sohoinst + '_diff'
    if "Multi" in mode: sohoinst = sohoinst + '_MGN'
    imagedir = os.path.join(imagedir, sohoinst)

    return sohoinst, imagedir