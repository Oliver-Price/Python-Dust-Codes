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
