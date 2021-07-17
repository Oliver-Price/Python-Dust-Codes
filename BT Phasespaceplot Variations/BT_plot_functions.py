#*************************************
#Various PIL methods for the BT module 
#*************************************
import numpy as np
from conversion_routines import round_to_base

#%%
#*************************
#DIRECT DUSTPLOT FUNCTIONS
#*************************

def r2xpix(r, border, pixwidth, rl, scale):
    return border*1.5 + (r - rl)*scale

def simt2xpix(simt, border, pixwidth, simtl, scale):
    return pixwidth + border*1.5 + (simtl - simt)*scale
    
def beta2ypix(beta, border, pixheight, betal, scale):
    return pixheight + border*1.5 + (betal - beta)*scale
    
def logbeta2ypix(beta, border, pixheight, betal, scale):
    return pixheight + border*1.5 + (np.log(betal) - np.log(beta))*scale

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

def logplotpixel(d,ta,ba,simres,dec,border,pixwt,pixhi,betal,simtl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = logbeta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
    t1 = simt2xpix(simres[ta,ba,0], border, pixwt, simtl, wscle)
    b2 = logbeta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
    t2 = simt2xpix(simres[ta,ba+1,0], border, pixwt, simtl, wscle)
    b3 = logbeta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
    t3 = simt2xpix(simres[ta+1,ba+1,0], border, pixwt,simtl, wscle)
    b4 = logbeta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
    t4 = simt2xpix(simres[ta+1,ba,0], border, pixwt, simtl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))

def logplotpixelr(d,ta,ba,simres,rvals,dec,border,pixwt,pixhi,betal,rl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = logbeta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
    t1 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b2 = logbeta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
    t2 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b3 = logbeta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
    t3 = r2xpix(rvals[ta+1], border, pixwt,rl, wscle)
    b4 = logbeta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
    t4 = r2xpix(rvals[ta+1], border, pixwt, rl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))

def plotpixelr(d,ta,ba,simres,rvals,dec,border,pixwt,pixhi,betal,rl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = beta2ypix(simres[ta,ba,1], border, pixhi, betal, hscle)
    t1 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b2 = beta2ypix(simres[ta,ba+1,1], border, pixhi, betal, hscle)
    t2 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b3 = beta2ypix(simres[ta+1,ba+1,1], border, pixhi, betal, hscle)
    t3 = r2xpix(rvals[ta+1], border, pixwt,rl, wscle)
    b4 = beta2ypix(simres[ta+1,ba,1], border, pixhi, betal, hscle)
    t4 = r2xpix(rvals[ta+1], border, pixwt, rl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))
    
def plotpixelbstart(d,ta,ba,simres,bstartvals,dec,border,pixwt,pixhi,betal,simtl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = beta2ypix(bstartvals[ta,ba], border, pixhi, betal, hscle)
    t1 = simt2xpix(simres[ta,ba,0], border, pixwt, simtl, wscle)
    b2 = beta2ypix(bstartvals[ta,ba+1], border, pixhi, betal, hscle)
    t2 = simt2xpix(simres[ta,ba+1,0], border, pixwt, simtl, wscle)
    b3 = beta2ypix(bstartvals[ta+1,ba+1], border, pixhi, betal, hscle)
    t3 = simt2xpix(simres[ta+1,ba+1,0], border, pixwt,simtl, wscle)
    b4 = beta2ypix(bstartvals[ta+1,ba], border, pixhi, betal, hscle)
    t4 = simt2xpix(simres[ta+1,ba,0], border, pixwt, simtl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))

def logplotpixelbstart(d,ta,ba,simres,bstartvals,dec,border,pixwt,pixhi,betal,simtl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = logbeta2ypix(bstartvals[ta,ba], border, pixhi, betal, hscle)
    t1 = simt2xpix(simres[ta,ba,0], border, pixwt, simtl, wscle)
    b2 = logbeta2ypix(bstartvals[ta,ba+1], border, pixhi, betal, hscle)
    t2 = simt2xpix(simres[ta,ba+1,0], border, pixwt, simtl, wscle)
    b3 = logbeta2ypix(bstartvals[ta+1,ba+1], border, pixhi, betal, hscle)
    t3 = simt2xpix(simres[ta+1,ba+1,0], border, pixwt,simtl, wscle)
    b4 = logbeta2ypix(bstartvals[ta+1,ba], border, pixhi, betal, hscle)
    t4 = simt2xpix(simres[ta+1,ba,0], border, pixwt, simtl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))
 
def logplotpixelrbstart(d,ta,ba,simres,bstartvals,rvals,dec,border,pixwt,pixhi,betal,rl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = logbeta2ypix(bstartvals[ta,ba], border, pixhi, betal, hscle)
    t1 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b2 = logbeta2ypix(bstartvals[ta,ba+1], border, pixhi, betal, hscle)
    t2 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b3 = logbeta2ypix(bstartvals[ta+1,ba+1], border, pixhi, betal, hscle)
    t3 = r2xpix(rvals[ta+1], border, pixwt,rl, wscle)
    b4 = logbeta2ypix(bstartvals[ta+1,ba], border, pixhi, betal, hscle)
    t4 = r2xpix(rvals[ta+1], border, pixwt, rl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))

def plotpixelrbstart(d,ta,ba,simres,bstartvals,rvals,dec,border,pixwt,pixhi,betal,rl,wscle,hscle,
                fillco1,fillco2,fillco3):
    b1 = beta2ypix(bstartvals[ta,ba], border, pixhi, betal, hscle)
    t1 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b2 = beta2ypix(bstartvals[ta,ba+1], border, pixhi, betal, hscle)
    t2 = r2xpix(rvals[ta], border, pixwt, rl, wscle)
    b3 = beta2ypix(bstartvals[ta+1,ba+1], border, pixhi, betal, hscle)
    t3 = r2xpix(rvals[ta+1], border, pixwt,rl, wscle)
    b4 = beta2ypix(bstartvals[ta+1,ba], border, pixhi, betal, hscle)
    t4 = r2xpix(rvals[ta+1], border, pixwt, rl, wscle)
    d.polygon([(t1,b1),(t2,b2),(t3,b3),(t4,b4)],fill=(fillco1,fillco2,fillco3,255))

def bt_setaxisup(simtu,simtl,betau,betal,logaxis,border,pixhi,hscle,pixwt,wscle):
    
    #Major divs in time
    tdivmajors = np.array([1.,2.,3.,4.,5.])
    tdivnos = (1/tdivmajors) * (simtu - simtl)
    tnodivs = 6
    tdividx = np.where((tdivnos <= tnodivs)==True)[0][0]
    tmajticks = np.arange(simtu, simtl-0.0001, -tdivmajors[tdividx])
    
    #Minor divs in time
    tdivminors = np.array([0.5,1,1,1,1])
    tminticks = np.arange(simtu, simtl-0.0001, -tdivminors[tdividx])
    tminticks = np.setdiff1d(tminticks,tmajticks)
      
    #Major divs in beta     
    bdivmajors = np.array([0.1,0.2,0.5,1,2,5,10,20,50])
    bdivnos = (1/bdivmajors) * (betau - betal)
    bnodivs = 10
    bdividx = np.where((bdivnos <= bnodivs)==True)[0][0]
    bu2majdv = round_to_base(betau, bdivmajors[bdividx])
    bl2majdv = round_to_base(betal, bdivmajors[bdividx])   
    bmajticks = np.arange(bu2majdv, bl2majdv-0.0001, -bdivmajors[bdividx])
    
    #Minor divs in beta  
    bdivminors = np.array([0.02,0.05,0.1,0.2,0.5,1,2,5,10])
    bu2mindv = round_to_base(betau, bdivminors[bdividx])
    bl2mindv = round_to_base(betal, bdivminors[bdividx])   
    bminticks = np.arange(bu2mindv, bl2mindv-0.0001, -bdivminors[bdividx])
    bminticks = np.setdiff1d(bminticks,bmajticks)
    
    if logaxis == True:
            bminlocs = logbeta2ypix(bminticks, border, pixhi, betal, hscle)
            bmajlocs = logbeta2ypix(bmajticks, border, pixhi, betal, hscle)
    else:
            bminlocs = beta2ypix(bminticks, border, pixhi, betal, hscle)
            bmajlocs = beta2ypix(bmajticks, border, pixhi, betal, hscle)
    
    tminlocs = simt2xpix(tminticks, border, pixwt, simtl, wscle)
    tmajlocs = simt2xpix(tmajticks, border, pixwt, simtl, wscle)

    return (bminlocs,bmajlocs,tminlocs,tmajlocs,tmajticks,bmajticks,tminticks,bminticks)

def rbr_setaxisup(ru,rl,betau,betal,logaxis,border,pixhi,hscle,pixwt,wscle):
    
    #Major divs in time
    tdivmajors = np.array([0.01,0.02,0.05,0.1,0.2])
    tdivnos = (1/tdivmajors) * (ru - rl)
    tnodivs = 6
    tdividx = np.where((tdivnos <= tnodivs)==True)[0][0]
    tu2majdv = round_to_base(ru, tdivmajors[tdividx])
    tl2majdv = round_to_base(rl, tdivmajors[tdividx])   
    tmajticks = np.arange(tu2majdv, tl2majdv-0.0001, -tdivmajors[tdividx])
    tmajticks = np.round(tmajticks,2)
    
    #Minor divs in time
    tdivminors = np.array([0.002,0.005,0.01,0.02,0.05])
    tu2mindv = round_to_base(ru, tdivminors[tdividx])
    tl2mindv = round_to_base(rl, tdivminors[tdividx])   
    tminticks = np.arange(tu2mindv, tl2mindv-0.0001, -tdivminors[tdividx])
    tminticks = np.round(tminticks,2)
    tminticks = np.setdiff1d(tminticks,tmajticks)
      
    #Major divs in beta     
    bdivmajors = np.array([0.1,0.2,0.5,1,2,5,10,20,50])
    bdivnos = (1/bdivmajors) * (betau - betal)
    bnodivs = 10
    bdividx = np.where((bdivnos <= bnodivs)==True)[0][0]
    bu2majdv = round_to_base(betau, bdivmajors[bdividx])
    bl2majdv = round_to_base(betal, bdivmajors[bdividx])   
    bmajticks = np.arange(bu2majdv, bl2majdv-0.0001, -bdivmajors[bdividx])
    bmajticks = np.clip(bmajticks,betal,betau)
    if (0 in bmajticks) and (logaxis == True):
        bmajticks[bmajticks == 0] = float('%.1g' % betal)
    
    #Minor divs in beta  
    bdivminors = np.array([0.02,0.05,0.1,0.2,0.5,1,2,5,10])
    bu2mindv = round_to_base(betau, bdivminors[bdividx])
    bl2mindv = round_to_base(betal, bdivminors[bdividx])   
    bminticks = np.arange(bu2mindv, bl2mindv-0.0001, -bdivminors[bdividx])
    bminticks = np.setdiff1d(bminticks,bmajticks)
    
    if logaxis == True:
            bminlocs = logbeta2ypix(bminticks, border, pixhi, betal, hscle)
            bmajlocs = logbeta2ypix(bmajticks, border, pixhi, betal, hscle)
    else:
            bminlocs = beta2ypix(bminticks, border, pixhi, betal, hscle)
            bmajlocs = beta2ypix(bmajticks, border, pixhi, betal, hscle)
    
    tminlocs = r2xpix(tminticks, border, pixwt, rl, wscle)
    tmajlocs = r2xpix(tmajticks, border, pixwt, rl, wscle)

    return (bminlocs,bmajlocs,tminlocs,tmajlocs,tmajticks,bmajticks,tminticks,bminticks)
        
#%% Antiquated
    
#scaling function to increase contrast
def greyscale_remap(upper_saturation,lower_saturation, mode = 'lin'):
    
    newmap = np.zeros((256),dtype=int)
    for col in range(upper_saturation, 256):
        newmap[col] = 255
    for col in range(0, lower_saturation):
        newmap[col] = 0
      
    if mode == 'lin':
        grad = 255/(upper_saturation - lower_saturation)   
        for col in range(lower_saturation, upper_saturation):
            newmap[col] = int(round((col - lower_saturation)*grad))
            
    if mode == 'log':
        print ('Using Log Scale')
        newmap[lower_saturation:upper_saturation+1] = np.rint(np.logspace( 
        np.log10(1), np.log10(255),
        num = (upper_saturation - lower_saturation + 1) ))
        
    return newmap
