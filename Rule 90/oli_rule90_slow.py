# -*- coding: utf-8 -*-

#import
import numpy as np
from scipy.ndimage.filters import convolve
import matplotlib.pyplot as plt
import time

def getbelow(left,middle,right):
    if left == 1:
        if middle == 1:
            if right == 1:
                below = 0
            elif right == 0:
                below = 1
        elif middle == 0:
            if right == 1:
                below = 0
            elif right == 0:
                below = 1
    elif left == 0:
        if middle == 1:
            if right == 1:
                below = 1
            elif right == 0:
                below = 0
        elif middle == 0:
            if right == 1:
                below = 1
            elif right == 0:
                below = 0  
    return below

#%%
    
#timing
start_time = time.clock() 

#generate
width = 100 ; timesteps = 100
array = np.zeros((width,timesteps),dtype=int)
array[width//2,0] = 1


for t in range(timesteps-1):
    left = 0
    middle = array[0,t]
    right = array[1,t]
    array[0,t+1] = getbelow(left,middle,right)
    x = 1
    while (x < width - 1):
       array[x,t+1] = getbelow(array[x-1,t],array[x,t],array[x+1,t]) 
       x += 1
    left = array[x-1,t]
    middle = array[x,t]
    right = 0
    array[x,t+1] = getbelow(left,middle,right)

#get time
print("--- %s seconds ---" % (time.clock() - start_time)) 

#plot
plt.imshow(array, aspect='auto', interpolation='none')