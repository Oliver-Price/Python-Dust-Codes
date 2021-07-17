#auxiliary functions for easy IO
import easygui
import datetime
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
    
        if not os.path.exists(os.path.join(imagedir,'temporary')):
            os.makedirs(os.path.join(imagedir,'temporary'))
    
        #making a 1d image for doing ra/dec coordinates
        fitstemp = os.path.join(os.path.join(imagedir,'temporary'), 'temporary_' + fitsinfile)
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
                 (('SOHO' in obslocstr) or ('Soho' in obslocstr)),
                 (('ISS' in obslocstr) or ('iss' in obslocstr))])
    name_locs = np.array(['Earth', 'Stereo_A', 'Stereo_B', 'Soho', 'ISS'])
    case_locs = np.size(np.nonzero(bool_locs))
    if case_locs > 1:
        obsmsg = "Please select observer location"
        obschoices = name_locs[bool_locs].tolist()
        obsloc = easygui.buttonbox(obsmsg, choices=obschoices)
        imagedir = os.path.join(imagedir, obsloc)
    elif case_locs == 1:
        obsloc = name_locs[bool_locs][0]
        imagedir = os.path.join(imagedir, obsloc)
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

    sname_locs = np.array(['Clear', 'Blue','Orange'])
    sohomsg = "Please select SOHO Filter"
    sohochoices = sname_locs.tolist()
    sohoinst = easygui.buttonbox(sohomsg, choices=sohochoices)
    sohoinst = 'C3_' + sohoinst
    
    mname_locs = np.array(['Base', 'Difference', 'Multiscale Gaussian Normalised'])
    modemsg = "Please select image type"
    modechoices = mname_locs.tolist()
    mode = easygui.buttonbox(modemsg, choices=modechoices)
    
    if "Diff" in mode: sohoinst = sohoinst + '_diff'
    if "Multi" in mode: sohoinst = sohoinst + '_MGN'
    imagedir = os.path.join(imagedir, sohoinst)

    return sohoinst, imagedir
    
def get_hih_low(comdenom,obsloc,inst):
    
    if comdenom == 'c2011l4':
        if obsloc == 'Stereo_B':
            low = 30000; hih = 800000
            if 'MGN' in inst:
                low = -0.25; hih = 0.5
            if 'diff' in inst:
                low = -270; hih = 240
        elif obsloc == 'Stereo_A':
            low = 3000; hih = 70000
        elif obsloc == 'Earth':
            low = 0; hih = 255
            
    elif comdenom == 'c2006p1':
        if "Stereo" in obsloc:
            low = 10000; hih = 1500000
            if 'diff' in inst: 
                low = -100; hih = 100
            elif 'MGN' in inst:  
                low = -0.7; hih = 1.25
        elif obsloc == 'Earth':
            low = 0; hih = 255
        elif obsloc == 'Soho':
            if 'Clear' in inst:
                low = 4.8e-10; hih = 1.4e-9
                if 'MGN' in inst:  
                    low = -0.6; hih = 0.8
            if  'Blue' in inst:
                low = 200; hih =  3500
                #low = 350; hih = 1200
                if 'MGN' in inst:  
                    low = -0.5; hih = 1.2
            if 'Orange' in inst:
                low = 300; hih = 3000
    
    elif comdenom == 'c2002v1':
        if obsloc == 'Earth':
            low = 0; hih = 255
        elif obsloc == 'Soho':
            low = 1e-13; hih = 1e-11
            if 'MGN' in inst:  
                low = 0; hih = 1
            elif 'diff' in inst:
                low = -7e-13; hih = 7e-13
                
    elif comdenom == 'c2011w3':
        if obsloc == 'Earth' or obsloc == 'ISS':
            low = 0; hih = 255
        elif obsloc == 'Soho':
            low = 1000; hih = 10000
            if 'MGN' in inst:  
                low = -0.2; hih = 1
            elif 'diff' in inst:
                low = -1000; hih = 1000
            
    elif comdenom == 'c1965s1':   
        low = 0; hih = 255
    
    elif comdenom == '96P':
        low = -550; hih = 380
    
    else:
        low = 0; hih = 255    
                
    return hih, low


