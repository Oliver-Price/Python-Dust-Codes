#********************************
#Comet Orbital Diagnostic Methods
#********************************

import numpy as np
from conversion_routines import pos2radec

def plot_orbit(comobs,comveceq,obsveceq,axisdata,d,comcel,trajfill,ra_img_lower,
               ra_img_higher,border,pixwidth,rafmin,scale,pixheight,decmin):
    
    #find rough cell range of traj from observer data
    obsraloc = np.where((comobs[:,5] > ra_img_lower) \
    & (comobs[:,5] < ra_img_higher))[0]
    obsdecloc = np.where((comobs[:,6] > axisdata[8]) \
    & (comobs[:,6] < axisdata[9]))[0]
    trajrough = np.intersect1d(obsraloc,obsdecloc)
    
    #use this to calculate ra and dec of comet for a purposely oversized range
    vno = 0; vext = int(np.size(trajrough))
    vtraj = np.empty((np.size(trajrough)+2*vext-1,11),dtype = float)
    tcellmax = min(trajrough[-1] + vext, np.shape(comveceq)[0])
    for tcell in xrange(trajrough[0] - vext,tcellmax):
        vtemp = comveceq[tcell,6:9] - obsveceq[comcel,6:9]    
        ptemp = pos2radec(vtemp)
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
    vtrajcel = np.where(abs(vtraj[:,5] - comcel) < 1e-4)[0][0]
    ltcomcel = vtraj[vtrajcel,4]

    #convert to ra and dec, and plot
    vtraj[:,2] = ra2xpix(vtraj[:,0],border,pixwidth,rafmin,scale)
    vtraj[:,3] = dec2ypix(vtraj[:,1],border,pixheight,decmin,scale)
    for ta in xrange(0, (np.shape(vtraj)[0]-1)):
        d.line([(vtraj[ta,2],vtraj[ta,3]),(vtraj[ta+1,2],vtraj[ta+1,3])],\
        fill = trajfill)
        
    return ltcomcel, vtraj, vtrajcel