#RK4 METHOD for simulating dust particle trajectories
'''
INPUTS: BETA, TIME TO SIMULATE FOR, DT,
        INPUT POSITION AND VELOCITIES, POSITION OF EARTH AT FINISHING TIME

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

def part_sim(beta, simt, dt, pstart, efinp, cor): #simt in days, dt in minutes
    simt = int(simt*1440) + cor
    traj = pstart
    for t in np.arange(0,simt-30,dt):
        traj = rk4(traj, beta, dt) #traj updated to t+1
    t += dt    
    dr = np.linalg.norm(traj[0:3] - efinp)
    tret = t + 8.316746397269274*dr    
    while (tret < simt):
        traj = rk4(traj, beta, dt) #traj updated to t+1
        t += dt
        dr = np.linalg.norm(traj[0:3] - efinp)
        tret = t + 8.316746397269274*dr
    return (t,traj)
    
def diff(pstate, beta):
    diffed = np.empty((6),dtype = float)
    diffed[0:3] = 1.1574074074074073e-05*pstate[3:6]
    invr = 1/(np.linalg.norm(pstate[0:3]))
    acc = (-3.4249098175024633e-09)*(1-beta)*invr*invr*invr
    diffed[3:6] = acc*pstate[0:3] #3rd invr normalises the position vector
    return diffed
    
def rk4(pstate, beta, dt): 
    dt = int(dt*60) #DT INPUT IN MINUTES
    k1 = diff(pstate, beta)*dt
    k2 = diff(pstate+k1*0.5, beta)*dt
    k3 = diff(pstate+k2*0.5, beta)*dt
    k4 = diff(pstate+k3, beta)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666
    
    