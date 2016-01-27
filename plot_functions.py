#**********************************************
#Various methods for images with the PIL module
#**********************************************
import numpy as np
import math as m

#%%
#***************************************
#FINSON PROBSTEIN OVERLAY PLOT FUNCTIONS
#***************************************

def ra2xpix(ra, border, pixwidth, ramin, scale):
    return pixwidth + border*1.5 - (ra - ramin)*scale
    
def dec2ypix(dec, border, pixheight, decmin, scale):
    return pixheight + border*1.5 - (dec-decmin)*scale

def setaxisup(ramax,ramin,decmax,decmin,border,pixheight,pixwidth,scale):
    
    #Tables with various sets of possible major/minor division sizes
    divmajors = np.array([0.01,0.02,0.05,0.1,0.2,0.5,1,2,5])
    divminors = np.array([0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1])
    divrecips = np.array([100,50,20,10,5,2,1,0.5,0.2])
    divminrecips = np.array([500,200,100,50,20,10,5,2,1])
    
    #finds the extent of bordered region in ra/dec space
    ramaxb = (ramax + border*0.5/scale)
    raminb = (ramin - border*0.5/scale)
    decmaxb = (decmax + border*0.5/scale)
    decminb = (decmin - border*0.5/scale)
    
    #works out how many major ticks each set of divisions would generate
    radivnos = divrecips * (ramaxb - raminb + border/scale)
    decdivnos = divrecips * (decmaxb - decminb + border/scale)
    
    nodivs = 6 #chooses largest number of ticks on either axis
    radividx = (np.abs(radivnos-nodivs)).argmin()
    decdividx = (np.abs(decdivnos-nodivs)).argmin()
    dividx = max(radividx, decdividx) #sets aspect ratio equal

    #rounds up/down from lowest/highest bordered values to find lowest/highest ticks
    rahidiv = m.floor(ramaxb*divrecips[dividx])*divmajors[dividx]
    ralodiv = m.ceil(raminb*divrecips[dividx])*divmajors[dividx]
    dechidiv = m.floor(decmaxb*divrecips[dividx])*divmajors[dividx]
    declodiv = m.ceil(decminb*divrecips[dividx])*divmajors[dividx]
    
    #as before, but for minor ticks
    rahimindiv = m.floor(ramaxb*divminrecips[dividx])*divminors[dividx]
    ralomindiv = m.ceil(raminb*divminrecips[dividx])*divminors[dividx]
    dechimindiv = m.floor(decmaxb*divminrecips[dividx])*divminors[dividx]
    declomindiv = m.ceil(decminb*divminrecips[dividx])*divminors[dividx]
    
    #generates full list of major/minor ticks for RA/DEC
    a = np.arange(ralodiv, rahidiv+1e-10, divmajors[dividx])
    c = np.arange(ralomindiv, rahimindiv+1e-10, divminors[dividx])
    e = np.arange(declodiv, dechidiv+1e-10, divmajors[dividx])
    g = np.arange(declomindiv, dechimindiv+1e-10, divminors[dividx])
    
    #converts these to position on image
    b = ra2xpix(a,border,pixwidth,ramin,scale)
    d = ra2xpix(c,border,pixwidth,ramin,scale)
    f = dec2ypix(e,border,pixheight,decmin,scale)
    h = dec2ypix(g,border,pixheight,decmin,scale)
    
    return(a,b,d,e,f,h,raminb,ramaxb,decminb,decmaxb)
    #don't need the actual values of minor ticks so ignore c and g

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