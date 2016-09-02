#*************************************
#Various PIL methods for the BT module 
#*************************************
import numpy as np

#%%
#*************************
#DIRECT DUSTPLOT FUNCTIONS
#*************************

def logsimt2xpix(simt, border, simtl, scale):
    return border*1.5 + (np.log10(simt) - np.log10(simtl))*scale
    
def linsimt2xpix(simt, border, simtl, scale):
    return border*1.5 + (simt - simtl)*scale
    
def beta2ypix(beta, border, pixheight, betal, scale):
    return pixheight + border*1.5 - (np.log10(beta) - np.log10(betal))*scale

#FUNCTION TO SLIM DOWN AN ARRAY OF RA AND DEC VALUES USING ANOTHER SET OF RA
#AND DEC VALUES AS REFRENCE

def radec_slim(ra_data, dec_data, ra_ref, dec_ref):

    #this assumes that the data arrays of ra and dec have the same shape
    rashape0 = np.shape(ra_data)[0]
    rashape1 = np.shape(ra_data)[1]
    
    raimin = np.min(ra_ref[np.nonzero(ra_ref)])
    raimax = np.max(ra_ref[np.nonzero(ra_ref)])
    decimin = np.min(dec_ref[np.nonzero(dec_ref)])
    decimax = np.max(dec_ref[np.nonzero(dec_ref)])
    
    distmin = abs(ra_data - raimin) + abs(dec_data - decimin)
    minloc = np.where(distmin == np.min(distmin))
    distmax = abs(ra_data - raimax) + abs(dec_data - decimax)
    maxloc = np.where(distmax == np.min(distmax))
    
    up_index_0 = max(minloc[0][0], maxloc[0][0])
    down_index_0 = min(minloc[0][0], maxloc[0][0])
    up_index_1 = max(minloc[1][0], maxloc[1][0])
    down_index_1 = min(minloc[1][0], maxloc[1][0])
    
    if up_index_0 != (rashape0 - 1):
        up_index_0 += 1
    if down_index_0 != 0:
        down_index_0 -= 1
    if up_index_1 != (rashape1 - 1):
        up_index_1 += 1
    if down_index_1 != 0:
        down_index_1 -= 1
        
    ra_slim = ra_data[down_index_0:up_index_0,down_index_1:up_index_1]
    dec_slim = dec_data[down_index_0:up_index_0,down_index_1:up_index_1]

    return(ra_slim, dec_slim, down_index_0, down_index_1)
    
#scaling function to increase contrast
def greyscale_remap(upper_saturation,lower_saturation, mode = 'lin'):
    
    newmap = np.zeros((256),dtype=int)
    for col in xrange(upper_saturation, 256):
        newmap[col] = 255
    for col in xrange(0, lower_saturation):
        newmap[col] = 0
      
    if mode == 'lin':
        grad = 255/(upper_saturation - lower_saturation)   
        for col in xrange(lower_saturation, upper_saturation):
            newmap[col] = int(round((col - lower_saturation)*grad))
            
    if mode == 'log':
        print 'Using Log Scale'
        newmap[lower_saturation:upper_saturation+1] = np.rint(np.logspace( 
        np.log10(1), np.log10(255),
        num = (upper_saturation - lower_saturation + 1) ))
        
    return newmap
