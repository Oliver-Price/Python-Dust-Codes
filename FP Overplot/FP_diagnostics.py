#********************************
#Comet Orbital Diagnostic Methods
#********************************

import numpy as np
from PIL import ImageFont
from conversion_routines import pos2radec
from FP_plot_functions import ra2xpix, dec2ypix

#%%*****************************
#Plot "true" orbit path of comet
#*******************************
def plot_orbit(comobs,comveceq,obsveceq,axisdata,d,comcel,trajfill,ra_img_lower,
               ra_img_higher,border,pixwidth,rafmin,scale,pixheight,decmin,fnt,
               featur_fill, con_in_image,fixwrapsbool):
    
    #find rough cell range of traj from observer data
    obsraloc = np.where((comobs[:,5] > ra_img_lower) \
    & (comobs[:,5] < ra_img_higher))[0]
    obsdecloc = np.where((comobs[:,6] > axisdata[8]) \
    & (comobs[:,6] < axisdata[9]))[0]
    trajrough = np.intersect1d(obsraloc,obsdecloc)
    
    if trajrough.size != 0:    
        
        #use this to calculate ra and dec of comet for a purposely oversized range
        vno = 0; vext = int(np.size(trajrough))
        vtraj = np.empty((np.size(trajrough)+2*vext-1,11),dtype = float)
        tcellmax = min(trajrough[-1] + vext, np.shape(comveceq)[0])
        for tcell in range(trajrough[0] - vext,tcellmax):
            vtemp = comveceq[tcell,6:9] - obsveceq[comcel,6:9]    
            ptemp = pos2radec(vtemp,fixwrapsbool)
            vtraj[vno,0] = ptemp[0]
            vtraj[vno,1] = ptemp[1]
            vtraj[vno,4] = tcell
            vtraj[vno,5] = tcell + round(np.linalg.norm(vtemp)*8.316746397269274)
            vtraj[vno,6:11] = comveceq[tcell,1:6]
            vno +=1
        
        #use these ra and dec values to slim data down to within image borders
        trajrange = np.intersect1d(
                    np.intersect1d( np.where(vtraj[:,0] < ra_img_higher)[0],
                                    np.where(vtraj[:,0] > ra_img_lower)[0]),
                    np.intersect1d( np.where(vtraj[:,1] < axisdata[9])[0],
                                    np.where(vtraj[:,1] > axisdata[8])[0]))                  
        vtraj = vtraj[trajrange[0]:trajrange[-1],:]
        
        #find relevant cell in vtraj and comveceq accounting for LT
        if con_in_image == True:
            vtrajcel = np.where(abs(vtraj[:,5] - comcel) < 1e-4)[0][0]
            ltcomcel = vtraj[vtrajcel,4]
        else: ltcomcel = None; vtrajcel = None
        
        #convert to ra and dec, and plot
        vtraj[:,2] = ra2xpix(vtraj[:,0],border,pixwidth,rafmin,scale)
        vtraj[:,3] = dec2ypix(vtraj[:,1],border,pixheight,decmin,scale)
        for ta in range(0, (np.shape(vtraj)[0]-1)):
            d.line([(vtraj[ta,2],vtraj[ta,3]),(vtraj[ta+1,2],vtraj[ta+1,3])],\
            fill = trajfill)
        
        #plot the path  
        d.text((2*border + pixwidth + 30,border + 200), \
        "Orbital Path:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 230),
        (2*border + pixwidth + 170,border + 230)], fill = trajfill)
        
    else: ltcomcel = None; vtraj = None; vtrajcel = None
        
    return ltcomcel, vtraj, vtrajcel
    
#%%***********************************
#Plot points along orbit path of comet
#*************************************
def plot_orbit_points(d,vtraj,smallfnt,featur_fill,pixheight,pixwidth,border):

    #find locations to plot on orbit
    mnsiz = 5
    cwid = 2
    max_orbit_points = 6
    orbit_cells = np.where(vtraj[:,10] == 0)[0]
    if np.size(orbit_cells) > max_orbit_points:
        orbit_cells = np.intersect1d(np.where(vtraj[:,9]%3 == 0)[0],orbit_cells)
    if np.size(orbit_cells) > max_orbit_points:
        orbit_cells = np.intersect1d(np.where(vtraj[:,9]%6 == 0)[0],orbit_cells)
    if np.size(orbit_cells) > max_orbit_points:
        orbit_cells = np.intersect1d(np.where(vtraj[:,9]%12 == 0)[0],orbit_cells) 
    if np.size(orbit_cells) > max_orbit_points:
        orbit_cells = np.intersect1d(np.where(vtraj[:,9]%24 == 0)[0],orbit_cells)
#    if np.size(orbit_cells) > max_orbit_points:
#        orbit_cells = np.intersect1d(np.where(vtraj[:,8]%2 == 0)[0],orbit_cells)   
#    if np.size(orbit_cells) > max_orbit_points:
#        orbit_cells = np.intersect1d(np.where(vtraj[:,8]%4 == 0)[0],orbit_cells)
#        
    #plot these locations on the orbit                                
    for midn in orbit_cells.tolist():
        d.line( [ (vtraj[midn,2] + mnsiz ,vtraj[midn,3]+ mnsiz ), 
        (vtraj[midn,2]- mnsiz ,vtraj[midn,3]- mnsiz ) ],
        fill = featur_fill , width= cwid)
        d.line([(vtraj[midn,2] - mnsiz ,vtraj[midn,3]+ mnsiz ),
        (vtraj[midn,2]+ mnsiz ,vtraj[midn,3]- mnsiz )],
        fill = featur_fill , width= cwid)
#        orbtxtx = 2; orbtxty = 2;
#        if (vtraj[midn,3] > 0.93*pixheight+1.5*border):
#            orbtxty = 0; orbtxtx = 3
#        if (vtraj[midn,2] > 0.93*pixwidth+1.5*border):
#            orbtxtx = -13
#        d.text((vtraj[midn,2] + orbtxtx*mnsiz,vtraj[midn,3] + orbtxty*mnsiz),
#        (str(int(vtraj[midn,6])) + '/' + str(int(vtraj[midn,7])) + '/' + 
#        str(int(vtraj[midn,8])) + "\n%02d:00") % int(vtraj[midn,9]) ,
#        font=smallfnt, fill= featur_fill)

#%%*******************
#Plot Sun-Earth Vector
#*********************
def plot_sunearth_vec(d,comveceq,obsveceq,ltcomcel,axisdata,ra_img_higher,
                      ra_img_lower,border,pixwidth,pixheight,rafmin,decmin,
                      scale,comsunfill,featur_fill,fnt,fixwrapsbool):
    
    #creates an array of ra/dec values along sun-comet line
    cmsam = 10001
    offsets = np.ones(cmsam) + np.linspace(-1,1,cmsam)
    comsun = np.empty((cmsam,4),dtype = float)
    for cs in range(0,cmsam):
        ctemp = pos2radec(offsets[cs]*comveceq[ltcomcel,6:9] - 
                            obsveceq[ltcomcel,6:9],fixwrapsbool)
        comsun[cs,0] = ctemp[0]
        comsun[cs,1] = ctemp[1]
    
    #slims this array down to within image limits
    csvrange = np.intersect1d(
               np.intersect1d( np.where(comsun[:,0] < ra_img_higher)[0],
                                np.where(comsun[:,0] > ra_img_lower)[0]),
               np.intersect1d( np.where(comsun[:,1] < axisdata[9])[0],
                                np.where(comsun[:,1] > axisdata[8])[0]))                  
    comsun = comsun[csvrange[0]:csvrange[-1],:]

    #convert to ra and dec, and plot
    comsun[:,2] = ra2xpix(comsun[:,0],border,pixwidth,rafmin,scale)
    comsun[:,3] = dec2ypix(comsun[:,1],border,pixheight,decmin,scale)
    for ca in range(0, (np.shape(comsun)[0]-1)):
        d.line([(comsun[ca,2],comsun[ca,3]),
                    (comsun[ca+1,2],comsun[ca+1,3])],
                    fill = comsunfill)
                    
    d.text((2*border + pixwidth + 30,border + 310), \
    "Comet-Sun\nVector:",font=fnt, fill= featur_fill)
    d.line([(2*border + pixwidth + 30,border + 355),
            (2*border + pixwidth + 170,border + 355)], fill = comsunfill) 

#%%****************************************************************
#Plot Uncertainty of comet's position along orbit (from imgtimehdr)
#******************************************************************
def plot_compos_unc(d,idls,vtraj,vtrajcel,border,pixwidth,fnt,trajucfill,featur_fill):

        sig_t = int(float(idls.sig_t))
        ura = vtraj [ (vtrajcel - sig_t) : (vtrajcel + sig_t) , 2 ]
        udec = vtraj [ (vtrajcel - sig_t) : (vtrajcel + sig_t) , 3 ]
        for ua in range(0, np.size(ura)-2):
            d.line([(ura[ua],udec[ua]),(ura[ua+1],udec[ua+1])],  
            fill = trajucfill)
            
        d.text((2*border + pixwidth + 30,border + 250), \
        "Orbital\nUncertainty:",font=fnt, fill= featur_fill)
        d.line([(2*border + pixwidth + 30,border + 294),
            (2*border + pixwidth + 170,border + 294)], fill = trajucfill)
 
#%%**********************
#Graph Annotation Methods        
#************************

def write_bt_ranges(d,border, pixwidth, fnt, featur_fill,
                    betau, betal, simtu, simtl, bt_anno_idx):

    d.text((2*border + pixwidth + 30,border + 420 + 50*bt_anno_idx), 
    ("Beta Range: " + '\n' + '  ' + str(betal) + ' to ' + str(betau)),
    font=fnt, fill= featur_fill)
    
    d.text((2*border + pixwidth + 30,border + 470 + 50*bt_anno_idx), 
    ("Ejection Time\nrange in days\nbefore image:" + '\n' + '  '
    + str(simtl) + ' to ' + str(simtu)),
    font=fnt, fill= featur_fill)

def write_properties(d,border, pixwidth, fnt, featur_fill, obsloc, filebase,
                     comcel, ctime, comobs):

    #display orbit plane angle
    cimgopa = comobs[comcel,7]
    d.text((2*border + pixwidth + 30,border), \
    "Plane Angle:" + '\n' + '  ' + "%.2f" % cimgopa ,
    font=fnt, fill= featur_fill)
    
    #display image author
    if obsloc == 'Earth':
        image_from = filebase.split('_')[-1]
    else:
        image_from = obsloc
    d.text((2*border + pixwidth + 30,border + 50), \
    "Image From: " + '\n' + '  ' + image_from, font=fnt, fill= featur_fill)
    
    #display date
    d.text((2*border + pixwidth + 30,border + 100), \
    "Image Date: " + '\n' + '  ' + ctime.isot[0:10] ,
    font=fnt, fill= featur_fill)

    #display time
    d.text((2*border + pixwidth + 30,border + 150), \
    "Image Time: " + '\n' + '  ' + ctime.isot[11:16] ,
    font=fnt, fill= featur_fill)    
