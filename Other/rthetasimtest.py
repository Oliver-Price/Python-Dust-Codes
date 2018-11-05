# -*- coding: utf-8 -*-
import numpy as np

def diff_rtheta(pstate):
    diff0 = pstate[2]
    diff1 = pstate[0]*pstate[3]
    diff3 = pstate[0]*pstate[3]*pstate[3]
    diff4 = 2*pstate[3]*pstate[2]
    return np.array([diff0,diff1,diff3,diff4])

def rk4_rtheta(pstate, dt):
    k1 = diff_rtheta(pstate)*dt
    k2 = diff_rtheta(pstate+k1*0.5)*dt
    k3 = diff_rtheta(pstate+k2*0.5)*dt
    k4 = diff_rtheta(pstate+k3)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666

def diff_xyz(pstate):
    diff0 = pstate[3]
    diff1 = pstate[4]
    diff2 = pstate[5]
    return np.array([diff0,diff1,diff2,0,0,0])   
    	
def rk4_xyz(pstate, dt):
    k1 = diff_xyz(pstate)*dt
    k2 = diff_xyz(pstate+k1*0.5)*dt
    k3 = diff_xyz(pstate+k2*0.5)*dt
    k4 = diff_xyz(pstate+k3)*dt
    return pstate + (k1+k2*2+k3*2+k4)*0.16666666666666666

def get_out(traj): 
    traj_out = np.zeros(6)
    traj_out[0] = traj[0]*np.cos(traj[1])
    traj_out[1] = traj[0]*np.sin(traj[1])
    traj_out[3] = (traj[2]*np.cos(traj[1]) - traj[0]*traj[3]*np.sin(traj[1]))
    traj_out[4] = (traj[2]*np.sin(traj[1]) + traj[0]*traj[3]*np.cos(traj[1]))
    return traj_out

traj = np.array([10,0,0,0.1])
trajxyz = get_out(traj)

t = 0
simt = 1
dt = 0.1

while (t + dt < simt): #go to up to 30 minutes before current
    traj = rk4_rtheta(traj, dt) #traj updated
    trajxyz = rk4_xyz(trajxyz, dt)
    t += dt
    
trajxyz_test = get_out(traj)