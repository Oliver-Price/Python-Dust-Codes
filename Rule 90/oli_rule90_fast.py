# -*- coding: utf-8 -*-

#import
import numpy as np
from scipy.ndimage.filters import convolve
import matplotlib.pyplot as plt
import time

#timing
start_time = time.clock() 

#generate
width = 100 ; timesteps = 100
array = np.zeros((width,timesteps),dtype=float)
array[width//2,0] = 1
kernel = [-1,0,1]
for t in range(timesteps-1):
    array[:,t+1] = abs(convolve(array[:,t],kernel, mode='constant', cval=0))

#get time
print("--- %s seconds ---" % (time.clock() - start_time)) 

#plot
plt.imshow(array.T, aspect='equal', interpolation='none')

#print("--- %s seconds ---" % (time.clock() - start_time)) 