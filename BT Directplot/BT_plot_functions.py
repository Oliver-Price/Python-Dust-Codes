#*************************************
#Various PIL methods for the BT module 
#*************************************
import numpy as np

#%%
#*************************
#DIRECT DUSTPLOT FUNCTIONS
#*************************
    
def simt2xpix(simt, border, pixwidth, simtl, scale):
    return pixwidth + border*1.5 + (simtl - simt)*scale
    
def beta2ypix(beta, border, pixheight, betal, scale):
    return pixheight + border*1.5 + (betal - beta)*scale

def plotpixel(d,ta,ba,simres,dec,border,pixwt,pixhi,betal,simtl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = beta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
    t1 = simt2xpix(simres[ta,ba,0], border, pixwt, simtl, wscle)
    b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
    t2 = simt2xpix(simres[ta,ba+1,0], border, pixwt, simtl, wscle)
    b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
    t3 = simt2xpix(simres[ta+1,ba+1,0], border, pixwt,simtl, wscle)
    b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
    t4 = simt2xpix(simres[ta+1,ba,0], border, pixwt, simtl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))


#%% Antiquated
    
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
