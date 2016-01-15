#************************************************************
#Various methods for plotting fits images with the PIL module
#************************************************************
import numpy as np
import math as m

def ra2ypix(ra, border, ramin, scale):
    return border*1.5 + (ra - ramin)*scale
    
def dec2xpix(dec, border, pixheight, decmin, scale):
    return pixheight + border*1.5 - (dec-decmin)*scale

def setaxisup(ramax,ramin,decmax,decmin,border,pixheight,scale):
    
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
    b = ra2ypix(a,border,ramin,scale)
    d = ra2ypix(c,border,ramin,scale)
    f = dec2xpix(e,border,pixheight,decmin,scale)
    h = dec2xpix(g,border,pixheight,decmin,scale)
    
    return(a,b,d,e,f,h,raminb,ramaxb,decminb,decmaxb)
    #don't need the actual values of minor ticks so ignore c and g
