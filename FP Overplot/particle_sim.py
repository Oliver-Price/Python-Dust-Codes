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
'''
import numpy as np
import math as m

def part_sim(beta, simt, ndt, ltdt, pstart, efinp, cor): 
    simt = simt*10 + cor #simt in minutes
    t = 0
    dt = (simt-30)/ndt
    traj = pstart
    for n in xrange(0,ndt,1): #go to up to 30 minutes before current
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

def frag_sim(beta, bfac, simt, ftime, ndt, ltdt, pstart, efinp, cor): 
    simt = simt*10 + cor #simt in minutes
    ftime = ftime*24*60
    t = 0
    dt = (4320)/ndt
    traj = pstart
    while (t <= ftime) and (t < simt-30):#go to up to 30 minutes before current
        traj = rk4(traj, beta, dt) #traj updated
        t += dt
    while (t > ftime) and (t < simt-30):
        traj = rk4(traj, beta*bfac, dt) #traj updated
        t += dt
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta*bfac, ltdt) #traj updated b the minute
        t += ltdt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)
    