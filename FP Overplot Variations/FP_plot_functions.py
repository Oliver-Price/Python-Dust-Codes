
#***************************************************
#Various PIL methods for the FINSON PROBSTEIN module 
#***************************************************
import numpy as np
import math as m
import matplotlib.pyplot as plt


#%%
#***********************
#FP Image plot functions
#***********************

def ra2xpix(ra, border, pixwidth, ramin, scale):
    return pixwidth + border*1.5 - (ra - ramin)*scale
    
def dec2ypix(dec, border, pixheight, decmin, scale):
    return pixheight + border*1.5 - (dec-decmin)*scale
    
def plotpixel(d,x,y,ra_m,dec,border,pixwidth,pixheight,decmin,rafmin,scale,
                fillco1,fillco2,fillco3):
    ra1 = ra2xpix(ra_m[x,y],border,pixwidth,rafmin,scale)
    ra2 = ra2xpix(ra_m[x+1,y],border,pixwidth,rafmin,scale)
    ra3 = ra2xpix(ra_m[x+1,y+1],border,pixwidth,rafmin,scale)
    ra4 = ra2xpix(ra_m[x,y+1],border,pixwidth,rafmin,scale)
    dec1 = dec2ypix(dec[x,y],border,pixheight,decmin,scale)
    dec2 = dec2ypix(dec[x+1,y],border,pixheight,decmin,scale)
    dec3 = dec2ypix(dec[x+1,y+1],border,pixheight,decmin,scale)
    dec4 = dec2ypix(dec[x,y+1],border,pixheight,decmin,scale)
    d.polygon([(ra1,dec1),(ra2,dec2),(ra3,dec3),(ra4,dec4)] ,\
    fill=(fillco1,fillco2,fillco3,255))
    
def plotpixel2(d,x,y,rp,dp,fillco1,fillco2,fillco3):
    d.polygon([(rp[x,y],dp[x,y]),(rp[x+1,y],dp[x+1,y]),
               (rp[x+1,y+1],dp[x+1,y+1]),(rp[x,y+1],dp[x,y+1])] ,\
    fill=(fillco1,fillco2,fillco3,255))
        

#%%
#*****************
#FP Axis functions
#*****************

#axisdata 0/1 - ra major divisions --- values and pixel locations
#axisdata 2 - ra minor divisions --- pixel locations
#axisdata 3/4 - dec major divisions --- values and pixel locations
#axisdata 5 - dec minor divisions --- pixel locations
#axisdata 6/7 min/max RA in border
#axisdata 8/9 min/max DEC in border

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
#**********************
#FP Grid Plot Functions
#**********************
    
def draw_syndynes(dynfill,d,simres,bno,rapixl,decpixl,tmin,tmax,bidx_list,spacing,bvals):
   
    print('Syndynes')
    for ba in np.append(bidx_list[:-1][::spacing],bidx_list[-1]).tolist():
        for ta in range(tmin[ba], tmax[ba]):
            d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
            (simres[ta+1,ba,12],simres[ta+1,ba,13])],\
            fill = dynfill)
        #print (bvals[bidx_list[ba]])
            
#        d.line([(rapixl,decpixl), \
#        (simres[0,ba,12],simres[0,ba,13])],\
#        fill = dynfill)

def draw_synchrones(chrfill,d,simres,tno,rapixl,decpixl,bmin,bmax,tidx_list,spacing,tvals,w = 1):

    print('Synchrones')
    for ta in np.append(tidx_list[:-1][::spacing],tidx_list[-1]).tolist():
        for ba in range(bmin[ta], bmax[ta]):
            d.line([(simres[ta,ba,12],simres[ta,ba,13]), \
            (simres[ta,ba+1,12],simres[ta,ba+1,13])],\
            fill = chrfill,width = w)
        #print (tvals[tidx_list[ta]])      
#        d.line([(rapixl,decpixl), \
#        (simres[ta,0,12],simres[ta,0,13])],\
#        fill = chrfill)            
            
def draw_datap(drfill,d,simres,xsiz = 2):

    all_points = np.where(simres[:,:,14]==1)
    for pidx in range(0,np.size(all_points[0])):
        ta = all_points[0][pidx]; ba = all_points[1][pidx]
        d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] - xsiz ) ,
                  ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] + xsiz ) ] ,
                  fill = drfill )  
        d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] + xsiz ) ,
                  ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] - xsiz ) ] ,
                  fill = drfill )

def draw_datap_CME(drfill,d,simres,xsiz = 2):

    all_points = np.where(simres[:,:,17]==1)
    for pidx in range(0,np.size(all_points[0])):
        ta = all_points[0][pidx]; ba = all_points[1][pidx]
        d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] - xsiz ) ,
                  ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] + xsiz ) ] ,
                  fill = drfill )  
        d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] + xsiz ) ,
                  ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] - xsiz ) ] ,
                  fill = drfill )
                      
def draw_data_reg(drfill,d,simres,bmax,bmin,bidx_list,tmax,tmin,tidx_list,border,pixwidth,lwidth = 5):                
    
    for ba in bidx_list[:-1].tolist():
        d.line( [ ( simres[tmax[ba],ba,12]  , simres[tmax[ba],ba,13] ) ,
                      ( simres[tmax[ba+1],ba+1,12] , simres[tmax[ba+1],ba+1,13] ) ]
                      , fill = drfill , width = lwidth)
        d.line( [ ( simres[tmin[ba],ba,12]  , simres[tmin[ba],ba,13] ) ,
                      ( simres[tmin[ba+1],ba+1,12] , simres[tmin[ba+1],ba+1,13] ) ]
                      , fill = drfill , width = lwidth)
    for ta in tidx_list[:-1].tolist():
        d.line( [ ( simres[ta,bmax[ta],12]  , simres[ta,bmax[ta],13] ) ,
                      ( simres[ta+1,bmax[ta+1],12] , simres[ta+1,bmax[ta+1],13] ) ]
                      , fill = drfill , width = lwidth)
        d.line( [ ( simres[ta,bmin[ta],12]  , simres[ta,bmin[ta],13] ) ,
                      ( simres[ta+1,bmin[ta+1],12] , simres[ta+1,bmin[ta+1],13] ) ]
                      , fill = drfill , width = lwidth)
                      
def draw_phase_points(d,simres,xsiz = 2):

    all_points = np.where(simres[:,:,14]==1)
    colormap = (plt.cm.winter(np.linspace(0,1,101))*255).astype(int)
    smin = np.min(simres[:,:,17][np.nonzero(simres[:,:,17])])
    smax = np.max(simres[:,:,17])
    grad = smax - smin
    
    for pidx in range(0,np.size(all_points[0])):
        ta = all_points[0][pidx]; ba = all_points[1][pidx]
        cval = int((simres[ta,ba,17]-smin)/grad*100)
        d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] - xsiz ) ,
                  ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] + xsiz ) ] ,
                  fill = tuple(colormap[cval]) )  
        d.line( [ ( simres[ta,ba,12] - xsiz , simres[ta,ba,13] + xsiz ) ,
                  ( simres[ta,ba,12] + xsiz , simres[ta,ba,13] - xsiz ) ] ,
                  fill = tuple(colormap[cval]) )
    
    return (smax, smin, grad, colormap)
#%%
#****************************
#FP Grid Annotation Functions
#****************************

def annotate_plotting(d,drawopts,border,pixwidth,fnt,featur_fill,dynfill,chrfill,drfill):

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
                
    elif (drawopts == "Data Region Enclosed"):
        bt_anno_idx = 0 
        d.text((2*border + pixwidth + 30,border + 370), \
                "Data Region:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 405),
                (2*border + pixwidth + 170,border + 405)], fill = drfill) 

    else: bt_anno_idx = 0
    return bt_anno_idx

def annotate_dustphase(d,border,pixwidth,featur_fill,fnt,smax,smin,grad,colormap):
   
    d.text((2*border + pixwidth + 30,border + 370), "Phase Angle:",font=fnt, fill= featur_fill)   
   
    for i in range(0,101):
        d.line([(2*border + pixwidth + 30,border + 405 + 2*i),
                (2*border + pixwidth + 100,border + 405 + 2*i)],
                fill = tuple(colormap[i]))
        d.line([(2*border + pixwidth + 30,border + 405 + 2*i + 1),
                (2*border + pixwidth + 100,border + 405 + 2*i + 1)],
                fill = tuple(colormap[i])) 
                
    for l in range(0,101,25):
        d.text((2*border + pixwidth + 110,border + 395 + 2*l), \
                ("%.0f" % (smin+grad*l/100)), font=fnt, fill= featur_fill)

    return 4

#%%
    
def getdustphase(sundusvec,sunobsvec):
    obsdusvec = sundusvec - sunobsvec
    dprod = obsdusvec[0]*sundusvec[0] + obsdusvec[1]*sundusvec[1] + obsdusvec[2]*sundusvec[2]
    obsdusmag = m.sqrt(obsdusvec[0]*obsdusvec[0] + obsdusvec[1]*obsdusvec[1] + obsdusvec[2]*obsdusvec[2])
    sundusmag = m.sqrt(sundusvec[0]*sundusvec[0] + sundusvec[1]*sundusvec[1] + sundusvec[2]*sundusvec[2])
    return m.acos(dprod/obsdusmag/sundusmag)*180/m.pi