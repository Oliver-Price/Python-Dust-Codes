#%%
import numpy as np
from scipy.ndimage.filters import convolve
import matplotlib.pyplot as plt
import time

getbelow_arr = np.zeros((2,2,2),dtype=int)
getbelow_arr[1,1,0] = 1
getbelow_arr[0,1,1] = 1
getbelow_arr[1,0,0] = 1
getbelow_arr[0,0,1] = 1

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
    array[0,t+1] = getbelow_arr[left,middle,right]
    x = 1
    while (x < width - 1):
       array[x,t+1] = getbelow_arr[array[x-1,t],array[x,t],array[x+1,t]]
       x += 1
    left = array[x-1,t]
    middle = array[x,t]
    right = 0
    array[x,t+1] = getbelow_arr[left,middle,right]

#get time
print("--- %s seconds ---" % (time.clock() - start_time)) 

#plot
plt.imshow(array.T, aspect='auto', interpolation='none')
