#RK4 METHOD for simulating dust particle trajectories
'''
INPUTS: BETA, TIME TO SIMULATE FOR, DT,
        INPUT POSITION AND VELOCITIES, POSITION OF EARTH AT FINISHING TIME

simt in 10s of minutes , ndt is now number of dt
ltdt is the dt when calculated lt correction

Diff function does this:
Input = Particle[r,v]	(Array with position and velocity)
dr = v;
a = -GM/|r|*|r|;
dv = a dot R(comet->sun) (Project a radially)
return dParticle = Array[dr,dv];

Distance is in AU with speed in AU/day but the simulation runs in minutes.

For the RK4 method we convert minutes to seconds then:
1.1574074074074073e-05 converts AU/day to AU/seconds

We also need solar GM in units of AU^3/day/second
Solar GM = 1.32712440018e20 m^3 s^-2
3.4249098175024633e-09 = 1.32712440018e20 / 149597870700^3 * 86400

For the output we assume the earth doesn't move a large distance away from
comet at end of it's traj. The final position is given where the simt is
reached accounting for a LT time to earth.
8.316746397269274 is c in Au/minute.

OUTPUTS: Final location and finishing time

NOTE:
    This method use a coarse 10 minute interval for the input for the position of the comet e.g. comvec10
    This means that comets with very extensive and old dust clouds (e.g. hale-bopp) can be measured
    However it's OTT for the majority of the cases, so I mainly use the part_sim_fine version below
    You can use the original version here but make sure to include the time correction (dtmin in main script)
    This method expects a simt in 10s of minutes (weird, I know)
'''
import numpy as np
import math as m

def part_sim(beta, simt, nperday, ltdt, pstart, efinp, cor): 
    simt = simt*10 + cor #simt in minutes
    dt1=1440/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    while (t + dt + 20 < simt): #go to up to 30 minutes before current
        traj = rk4(traj, beta, dt) #traj updated
        t += dt
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)
    
def diff(pstate, beta):
    diff0 = 1.1574074074074073e-05*pstate[3]
    diff1 = 1.1574074074074073e-05*pstate[4]
    diff2 = 1.1574074074074073e-05*pstate[5]
    invr = 1/(m.sqrt(pstate[0]**2 + pstate[1]**2 + pstate[2]**2))
    acc = (-3.4249098175024633e-09)*(1-beta)*invr*invr*invr
    diff3 = acc*pstate[0] #3rd invr normalises the position vector
    diff4 = acc*pstate[1]
    diff5 = acc*pstate[2] 
    return np.array([diff0,diff1,diff2,diff3,diff4,diff5])
    
def rk4(pstate, beta, dt): 
    dt = dt*60 #DT INPUT IN MINUTES
    k1 = diff(pstate, beta)*dt
    k2 = diff(pstate+k1*0.5, beta)*dt
    k3 = diff(pstate+k2*0.5, beta)*dt
    k4 = diff(pstate+k3, beta)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666


#%%
'''
This Method uses the normal comveceq position
'''

def part_sim_fine(beta, simt, nperday, ltdt, pstart, efinp): 
    dt1=1440/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    while (t + dt + 20 < simt): #go to up to 30 minutes before current
        traj = rk4(traj, beta, dt) #traj updated
        t += dt
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)

#%%tracking dust

def part_sim_fine_track(beta, simt, nperday, pstart): 
    dt1=1440/nperday
    dt2=float(simt)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    total = pstart
    times = np.zeros((1,1))
    while (t < simt):
        traj = rk4(traj, beta, dt) #traj updated
        t += dt
        total = np.vstack([total,traj])
        times = np.vstack([times,t])
    return (times,total)

#%%

'''
BASIC LORENTZ FORCE SIM

Same as before with lorentz force

Lorentz acc = (12*Eo/C^2)*(B*V*p*v/Q)*(Z*Ro/R^2)*Beta^2*(m/AU)*(day/s)

Constant = (12*Eo/C^2)*(m/AU)*(day/s)
With Ro as 1

Parameters = (B*V*p*v/Q)

Kramer:
B = 3e-9; V = 5, p = 1000, v = 5e5, Q = 1
parameters = 7.5

Lorentz acc = Constant * Parameters * (Z/R/R) * Beta^2

Constant of equation is 1.8495809331145113e-10
from 8.854187817e-12*12*86400/149597870700/(5.76e-4**2)

'''

def part_sim_lorentz(beta, simt, nperday, ltdt, pstart, efinp, cor, tswap, mag = 3e-9,  voltage = 2, density = 1000, swspeed =  3e5, qpr = 1, sign = 1): 
    pars = mag*voltage*density*swspeed/qpr
    tswap = tswap*1440
    simt = simt*10 + cor #simt in minutes
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    while (t + dt + 20 < simt) & (t < tswap): #go to up to 20 minutes before current
        #times.append(t)
        traj = rk4_lorentz(traj, beta, dt, pars,sign) #traj updated
        t += dt
    sign = -sign
    while (t + dt + 20 < simt): #go to up to 20 minutes before current
        traj = rk4_lorentz(traj, beta, dt, pars,sign) #traj updated
        t += dt
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr
    while (tret < simt):
        traj = rk4_lorentz(traj, beta, ltdt, pars,sign) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)
    
def diff_lorentz(pstate,beta,parameters,sign):
    diff0 = 1.1574074074074073e-05*pstate[3]
    diff1 = 1.1574074074074073e-05*pstate[4]
    diff2 = 1.1574074074074073e-05*pstate[5]
    invr = 1/(m.sqrt(pstate[0]**2 + pstate[1]**2 + pstate[2]**2))
    rho = m.sqrt(pstate[0]**2 + pstate[1]**2)
    gravrad = (-3.4249098175024633e-09)*(1-beta)*invr*invr*invr
    lorentz = sign*1.8495809331145113e-10*parameters*pstate[2]*invr*invr*invr*beta*beta
    diff3 = (gravrad + pstate[2]*lorentz/rho)*pstate[0]
    diff4 = (gravrad + pstate[2]*lorentz/rho)*pstate[1]
    diff5 = gravrad*pstate[2] - rho*lorentz
    print(pstate[2]*invr)
    return np.array([diff0,diff1,diff2,diff3,diff4,diff5])
    
def rk4_lorentz(pstate, beta, dt, parameters,sign): 
    dt = dt*60 #DT INPUT IN MINUTES
    k1 = diff_lorentz(pstate, beta,parameters,sign)*dt
    k2 = diff_lorentz(pstate+k1*0.5, beta,parameters,sign)*dt
    k3 = diff_lorentz(pstate+k2*0.5, beta,parameters,sign)*dt
    k4 = diff_lorentz(pstate+k3, beta,parameters,sign)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666

#%%

'''
NEAR SUN BETA BREAKDOWN SIM
'''

def part_sim_nearsun(beta, simt, nperday, ltdt, pstart, efinp, cor): 
    simt = simt*10 + cor #simt in minutes
    bdict = get_betadict()
    dt1=float(simt-30)*86400/simt/nperday
    dt2=float(simt-30)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    while (t + dt + 30 < simt): #go to up to 30 minutes before current
        traj = rk4_nearsun(traj, beta, bdict, dt) #traj updated
        t += dt
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4_nearsun(traj, beta, bdict, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)

def diff_nearsun(pstate, beta, bdict):
    diff0 = 1.1574074074074073e-05*pstate[3]
    diff1 = 1.1574074074074073e-05*pstate[4]
    diff2 = 1.1574074074074073e-05*pstate[5]
    r = m.sqrt(pstate[0]**2 + pstate[1]**2 + pstate[2]**2)
    invr = r**-1
    acc = (-3.4249098175024633e-09)*(1-beta*bdict.get(round(r,4),1))*invr*invr*invr
    diff3 = acc*pstate[0] #3rd invr normalises the position vector
    diff4 = acc*pstate[1]
    diff5 = acc*pstate[2] 
    return np.array([diff0,diff1,diff2,diff3,diff4,diff5])
    
def rk4_nearsun(pstate, beta, bdict, dt): 
    dt = dt*60 #DT INPUT IN MINUTES
    k1 = diff_nearsun(pstate, beta, bdict)*dt
    k2 = diff_nearsun(pstate+k1*0.5, beta, bdict)*dt
    k3 = diff_nearsun(pstate+k2*0.5, beta, bdict)*dt
    k4 = diff_nearsun(pstate+k3, beta, bdict)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666

def get_betadict():
    filename = r'C:\Users\Ollie\Documents\PhD\Comet_Data\Python Save Files\betafraction.dat'
    with open(filename, "r") as text_file:
        dat = text_file.readlines()
        a = np.zeros((len(dat)),dtype=float)
        b = np.zeros((len(dat)),dtype=float)
        for x in range(0,len(dat)):
            a[x] = float(dat[x][0:14].strip())
            b[x] = float(dat[x][14:-2].strip())
    return dict(zip(a,b))

#%%

'''
BASIC Sekanina-Farrell Striation Model

This simulation incorporates a change in beta by a factor of bfac after
a number of days equal to ftime
'''

def frag_sim_basic(beta, bsta, simt, ftime, nperday, ltdt, pstart, efinp, cor): 
    simt = simt*10 + cor #simt in minutes
    ftime = ftime*24*60 #frag time in days
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    while (t + dt + 20 < simt) and (t <= ftime):
        traj = rk4(traj, bsta, dt) #traj updated
        t += dt
    while (t + dt + 20 < simt) and (t > ftime): 
        traj = rk4(traj, beta, dt) #traj updated
        t += dt   
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)

'''
Fragmentation after a given exposure time
'''

def frag_sim_exposure(beta, bsta, simt, max_expos, nperday, ltdt, pstart, efinp, cor, tswap):
    simt = simt*10 + cor #simt in minutes
    tswap = tswap*1440
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    exposure_t = 0
    traj = pstart
    invr_old = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
    while (t + dt + 20 < simt) and (exposure_t <= max_expos):
        traj = rk4(traj, bsta, dt) #traj updated
        t += dt
        invr_new = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
        exposure_t += 0.00034722222222222224*dt*(invr_old*invr_old + invr_new*invr_new)
        invr_old = invr_new
    while (t + dt + 20 < simt): 
        traj = rk4(traj, beta, dt) #traj updated
        t += dt
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta, ltdt) #traj updated by the ltdt
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj,exposure_t)

#%%

'''
Including exposure time and lorentz force
'''

def frag_sim_exposure_lorentz(beta, bsta, simt, max_expos, nperday, ltdt, pstart, efinp, cor, tswap, mag = 3e-9, voltage = 2, density = 800, swspeed =  3e5, qpr = 1, sign = 1): 
    pars = mag*voltage*density*swspeed/qpr
    simt = simt*10 + cor #simt in minutes
    tswap = tswap*1440
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    exposure_t = 0
    traj = pstart
    invr_old = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
    while (t + dt + 20 < simt) and (exposure_t <= max_expos):
        traj = rk4_lorentz(traj, bsta, dt, pars,sign) #traj updated
        t += dt
        invr_new = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
        exposure_t += 0.00034722222222222224*dt*(invr_old*invr_old + invr_new*invr_new)
        invr_old = invr_new
        if t > tswap:
            sign = -sign
    while (t + dt + 20 < simt): 
        traj = rk4_lorentz(traj, beta, dt, pars,sign) #traj updated
        t += dt
        if t > tswap:
            sign = -sign
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4_lorentz(traj, beta, ltdt, pars,sign) #traj updated by the ltdt
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
        if t > tswap:
            sign = -sign
    return (t,traj,exposure_t)

#%%

'''
FINITE LIFETIME MODEL SIM w/lorentz force
'''

def F_s(S):
    return 0.48*10**(-1.3*(np.log10(S/0.24)**2))

def FLM_sim(bfrag, simt, nperday, ltdt, pstart, efinp, cor, tswap, h, c, s0, sc, mag = 3e-9, voltage = 2, density = 800, swspeed =  3e5, qpr = 1, sign = 1): 
   
    #initialisations
    pars = mag*voltage*density*swspeed/qpr
    simt = simt*10 + cor #simt in minutes
    tswap = tswap*1440
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)

    #parameters set up
    t = 0
    exposure_t = 0
    s = s0
    traj = pstart
    invr_old = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
    beta = h*F_s(s)
    
    while (t + dt + 20 < simt) and (s >= sc):
        
        #traj update
        traj = rk4_lorentz(traj, beta, dt, pars,sign) #traj updated
        t += dt
        
        #updating exposure time
        invr_new = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
        exposure_t += 0.00034722222222222224*dt*(invr_old*invr_old + invr_new*invr_new)
        invr_old = invr_new
        
        #updating s and beta
        s = s0*np.exp(-exposure_t/c)
        beta = h*F_s(s)
        
        if t > tswap:
            sign = -sign
            
    while (t + dt + 20 < simt): 
        
        traj = rk4_lorentz(traj, bfrag, dt, pars,sign) #traj updated
        t += dt
        
        #updating exposure time
        invr_new = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
        exposure_t += 0.00034722222222222224*dt*(invr_old*invr_old + invr_new*invr_new)
        invr_old = invr_new
        
        if t > tswap:
            sign = -sign
            
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    
    while (tret < simt):
        
        traj = rk4_lorentz(traj, bfrag, ltdt, pars,sign) #traj updated by the ltdt
        t += ltdt
        
        #updating exposure time
        invr_new = 1/(m.sqrt(traj[0]**2 + traj[1]**2 + traj[2]**2))
        exposure_t += 0.00034722222222222224*dt*(invr_old*invr_old + invr_new*invr_new)
        invr_old = invr_new
        
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
        
        if t > tswap:
            sign = -sign
            
    return (t,traj,exposure_t)

#%%

'''
Lorentz force diagnostic model for getting magnetic field diagnostics
'''

def part_sim_lorentz_diag(beta, simt, nperday, ltdt, pstart, efinp, cor, tswap, mag = 3e-9,  voltage = 2, density = 1000, swspeed =  3e5, qpr = 1, sign = 1): 
    a = []; times = []
    pars = mag*voltage*density*swspeed/qpr
    tswap = tswap*1440
    simt = simt*10 + cor #simt in minutes
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    t = 0
    traj = pstart
    while (t + dt + 20 < simt) & (t < tswap): #go to up to 20 minutes before current
        times.append(t)
        traj = rk4_lorentz(traj, beta, dt, pars,sign,a) #traj updated
        t += dt
    sign = -sign
    while (t + dt + 20 < simt): #go to up to 20 minutes before current
        times.append(t)
        traj = rk4_lorentz(traj, beta, dt, pars,sign,a) #traj updated
        t += dt
    print(t,dt,simt)
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr
    while (tret < simt):
        times.append(t)
        traj = rk4_lorentz(traj, beta, ltdt, pars,sign,a) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj,a,times)
    
def diff_lorentz_diag(pstate,beta,parameters,sign,a,b=False):
    diff0 = 1.1574074074074073e-05*pstate[3]
    diff1 = 1.1574074074074073e-05*pstate[4]
    diff2 = 1.1574074074074073e-05*pstate[5]
    invr = 1/(m.sqrt(pstate[0]**2 + pstate[1]**2 + pstate[2]**2))
    rho = m.sqrt(pstate[0]**2 + pstate[1]**2)
    gravrad = (-3.4249098175024633e-09)*(1-beta)*invr*invr*invr
    lorentz = sign*1.8495809331145113e-10*parameters*pstate[2]*invr*invr*invr*beta*beta
    diff3 = (gravrad + pstate[2]*lorentz/rho)*pstate[0]
    diff4 = (gravrad + pstate[2]*lorentz/rho)*pstate[1]
    diff5 = gravrad*pstate[2] - rho*lorentz
    if b == True:
        a.append(sign*3e-9*pstate[2]*invr*invr)
    return np.array([diff0,diff1,diff2,diff3,diff4,diff5])
    
def rk4_lorentz_diag(pstate, beta, dt, parameters,sign,a): 
    dt = dt*60 #DT INPUT IN MINUTES
    k1 = diff_lorentz(pstate, beta,parameters,sign,a,True)*dt
    k2 = diff_lorentz(pstate+k1*0.5, beta,parameters,sign,a)*dt
    k3 = diff_lorentz(pstate+k2*0.5, beta,parameters,sign,a)*dt
    k4 = diff_lorentz(pstate+k3, beta,parameters,sign,a)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666

#%%

'''
Lorentz force diagnostic model for getting magnetic field diagnostics
Note: Fairly sure this doesn't work
'''

def part_sim_rtheta(beta, simt, nperday, ltdt, pstart, efinp): 

	#set up timesteps
	dt1=float(simt-30)*86400/simt/nperday
	dt2=float(simt-30)/nperday    
	dt = min(dt1,dt2)

	#set up polar traj
	traj = np.empty(4)
	traj[0] = (pstart[0]**2 + pstart[1]**2)**0.5
	traj[1] = np.arctan2(pstart[1],pstart[0])%(2*np.pi)
	traj[2] = (pstart[0]*pstart[3] + pstart[1]*pstart[4])/traj[0]*1.1574074074074073e-05
	traj[3] = (-pstart[1]*pstart[3] + pstart[0]*pstart[4])/traj[0]/traj[0]*1.1574074074074073e-05

	#set up other things
	t= 0 
	mu = (1 - beta)
	
	#main simulation
	while (t + dt + 30 < simt): #go to up to 30 minutes before current
		traj = rk4_rtheta(traj, mu, dt) #traj updated
		t += dt
		  
	#up to light travel time
	dr = (efinp[0]**2 + traj[0]**2 + 2*efinp[0]*traj[0]*np.cos(efinp[1]-traj[1]))**0.5
	tret = t + 8.316746397269274*dr    
	while (tret < simt):
		traj = rk4_rtheta(traj, mu, ltdt) #traj updated b the minute
		t += ltdt
		dr = (efinp[0]**2 + traj[0]**2 + 2*efinp[0]*traj[0]*np.cos(efinp[1]-traj[1]))**0.5
		tret = t + 8.316746397269274*dr
	
	#back to xyz
	traj_out = np.zeros(6)
	traj_out[0] = traj[0]*np.cos(traj[1])
	traj_out[1] = traj[0]*np.sin(traj[1])
	traj_out[3] = (traj[2]*np.cos(traj[1]) - traj[0]*traj[3]*np.sin(traj[1]))*86400
	traj_out[4] = (traj[2]*np.sin(traj[1]) + traj[0]*traj[3]*np.sin(traj[1]))*86400
		
	return (t,traj_out)
    
def diff_rtheta(pstate, mu):
	diff0 = pstate[2]
	diff1 = pstate[0]*pstate[3]
	invr = 1/pstate[0]
	acc = (-3.9638468e-14)*mu*invr*invr
	diff3 = acc - pstate[0]*pstate[3]*pstate[3]
	diff4 = 2*pstate[3]*pstate[2]
	return np.array([diff0,diff1,diff3,diff4])   
	 
def rk4_rtheta(pstate, mu, dt): 
	dt = dt*60 #DT INPUT MINUTES TO SECONDS
	k1 = diff_rtheta(pstate, mu)*dt
	k2 = diff_rtheta(pstate+k1*0.5, mu)*dt
	k3 = diff_rtheta(pstate+k2*0.5, mu)*dt
	k4 = diff_rtheta(pstate+k3, mu)*dt
	return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666

#%%

'''
Lorentz force diagnostic model for getting magnetic field diagnostics
'''
import matplotlib.path as mplPath

def part_sim_CME(beta, simt, nperday, ltdt, pstart, efinp, cmedelta):
    
    dt1=float(simt-21)*86400/simt/nperday
    dt2=float(simt-21)/nperday
    dt = min(dt1,dt2)
    
    mult = {True:10,False:1}
    cmeSpeed = 888*0.000577548327
    		
    cmeMiddle = np.radians(310) - np.pi
    cmeEdge1 = cmeMiddle + np.radians(30)
    cmeEdge2 = cmeMiddle - np.radians(30)
    
    #rs = []
    #thts = []
    #crs = []
    
    t = 0    
    traj = pstart
    x = 0 

    while (t + dt + 20 < simt): #go to up to 30 minutes before current
        
        if x == 1: print(traj)
        
        cmeR = (cmedelta + t/1440 - simt/1440)*cmeSpeed
        cmeBack = cmeR - 0.03
        
        CME_box =((cmeR,cmeEdge1),
                  (cmeR,cmeEdge2),
                  (cmeBack,cmeEdge2),
                  (cmeBack,cmeEdge1),
                  (cmeR,cmeEdge1))
        
        com_Path = mplPath.Path(CME_box)  
        
        rDust = (traj[0]**2 + traj[1]**2 + traj[2]**2)**0.5
        thtDust = np.arctan2(traj[0],traj[1])%(2*np.pi)
        
        #crs.append(cmeR)
        #rs.append(rDust)
        #thts.append(thtDust)
        
        inCME = com_Path.contains_point((rDust,thtDust))

        betaUse = beta*mult[inCME]
        
        traj = rk4(traj, betaUse, dt) #traj updated
        t += dt
        
        x += 1
        
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    
    while (tret < simt):
        
        
        cmeR = (cmedelta + t/1440- simt/1440)*cmeSpeed
        cmeBack = cmeR - 0.03
        
        CME_box =((cmeR,cmeEdge1),
                  (cmeR,cmeEdge2),
                  (cmeBack,cmeEdge2),
                  (cmeBack,cmeEdge1),
                  (cmeR,cmeEdge1))
        
        com_Path = mplPath.Path(CME_box)  
        
        rDust = (traj[0]**2 + traj[1]**2 + traj[2]**2)**0.5
        thtDust = np.arctan2(traj[0],traj[1])%(2*np.pi)
        
        inCME = com_Path.contains_point((rDust,thtDust))
        
        betaUse = beta*mult[inCME]
        
        traj = rk4(traj, betaUse, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
        
    return (t,traj)
