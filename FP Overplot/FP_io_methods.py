#auxiliary functions for easy IO
import easygui
from astropy.io import fits
import os
import numpy as np
import sys
 
#%% - methods deals with various size data cubes
  
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