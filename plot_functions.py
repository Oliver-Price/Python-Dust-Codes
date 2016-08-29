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
    divmajors = np.array([0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20])
    divminors = np.array([0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5])
    divrecips = np.array([100,50,20,10,5,2,1,0.5,0.2,0.1,0.05])
    divminrecips = np.array([500,200,100,50,20,10,5,2,1,0.5,0.2])
    
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

#%%
#**********************
#FP Grid Plot Functions
#**********************
    
def draw_synchrones(dynfill,d,simres,bno,rapixl,decpixl,tmax):
    
    for ba in xrange(0, bno):
        d.line([(rapixl,decpixl), \
        (simres[0,ba,12],simres[0,ba,13])],\
        fill = dynfill)
        print 'BETA = ' + str(simres[0,ba,1])
        for ta in xrange(0, tmax[ba]):
            d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
            (simres[ta+1,ba,12],simres[ta+1,ba,13])],\
            fill = dynfill)

def draw_sydynes(chrfill,d,simres,tno,rapixl,decpixl,bmax):

    for ta in xrange(0, tno):
        d.line([(rapixl,decpixl), \
        (simres[ta,0,12],simres[ta,0,13])],\
        fill = chrfill)
        print 'TIME = ' + str(simres[ta,0,0])
        for ba in xrange(0, bmax[ta]):
            d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
            (simres[ta,ba+1,12],simres[ta,ba+1,13])],\
            fill = chrfill)
            
def draw_datap(drfill,d,simres,tno,bmax,xsiz = 2):

    for ta in xrange(0,tno):
        for ba in xrange(0, bmax[ta]):
            d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] - xsiz ) ,
                      ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] + xsiz ) ] ,
                      fill = drfill )  
            d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] + xsiz ) ,
                      ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] - xsiz ) ] ,
                      fill = drfill )
                      
def draw_data_reg(drfill,featur_fill,d,fnt,simres,border,pixwidth,width = 5):                
    
    b = d.line( [ ( simres[0,0,12]  , simres[0,0,13] ) ,
                  ( simres[0,bmax[0],12] , simres[0,bmax[0],13] ) ] ,
                    fill = drfill , width)
    b = d.line( [ ( simres[tno-1,0,12]  , simres[tno-1,0,13] ) ,
                  ( simres[tno-1,bmax[tno-1],12] ,
                   simres[tno-1,bmax[tno-1],13] ) ] ,    
                    fill = drfill , width)    
                    
    for ta in xrange(0,tno - 1):
        
        b = d.line( [ ( simres[ta,bmax[ta],12]  , simres[ta,bmax[ta],13] ) ,
                      ( simres[ta+1,bmax[ta+1],12] , simres[ta+1,bmax[ta+1],13] ) ]
                      , fill = drfill , width)
        b = d.line( [ ( simres[ta,0,12]  , simres[ta,0,13] ) ,
                      ( simres[ta+1,0,12] , simres[ta+1,0,13] ) ]
                      , fill = drfill , width)
                      
    d.text((2*border + pixwidth + 30,border + 370), \
            "Data Region:",font=fnt, fill= featur_fill)
    d.line([(2*border + pixwidth + 30,border + 405),
            (2*border + pixwidth + 170,border + 405)], fill = drfill)                  


#%%
#****************************
#FP Grid Annotation Functions
#****************************

def annotate_plotting(d,drawopts,border,pixwidth):

    if (drawopts == "Synchrones Only"):
        bt_anno_idx = 0
            
        d.text((2*border + pixwidth + 30,border + 370), \
                "Synchrones:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 405),
                (2*border + pixwidth + 170,border + 405)], fill = chrfill) 
            
    elif (drawopts == "Syndynes only"):
        bt_anno_idx = 0
        
        d.text((2*border + pixwidth + 30,border + 370), \
                "Syndynes:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 405),
                (2*border + pixwidth + 170,border + 405)], fill = dynfill) 
        
    elif (drawopts == "Synchrones and Syndynes"):
        bt_anno_idx = 1          
        
        d.text((2*border + pixwidth + 30,border + 370), \
                "Syndynes:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 405),
                (2*border + pixwidth + 170,border + 405)], fill = dynfill) 
                
        d.text((2*border + pixwidth + 30,border + 420), \
                "Synchrones:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 455),
                (2*border + pixwidth + 170,border + 455)], fill = chrfill) 
        
    elif (drawopts == "Synchrones, Syndynes and Data Points"):
        bt_anno_idx = 2     
        
        d.text((2*border + pixwidth + 30,border + 370), \
                "Syndynes:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 405),
                (2*border + pixwidth + 170,border + 405)], fill = dynfill) 
                
        d.text((2*border + pixwidth + 30,border + 420), \
                "Synchrones:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 455),
                (2*border + pixwidth + 170,border + 455)], fill = chrfill) 
                
        d.text((2*border + pixwidth + 30,border + 470), \
                "Data Points:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 505),
                (2*border + pixwidth + 170,border + 505)], fill = drfill)
                
    else: pass 

    return bt_anno_idx