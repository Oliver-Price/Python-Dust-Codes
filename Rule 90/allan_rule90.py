def shift(seq, n): # define a function to shift array either side for the xor on the original configuration
	n = n % len(seq)
	return seq[n:] + seq[:n]
import numpy as np
import random
import time
import matplotlib.pyplot as plt
start_time = time.clock() #time the whole process
x = [True]*100 # set up a list of booleans
x[50:] = [False]*50 #make half of them false

#for i,val in enumerate(x): # iterate over my list
#  x[i] = random.choice([True,False]) #make each entry random to start

plt.plot(np.array(list(map(int, x)))*2.5-9,'x')# plot the initial value
for i in range(0,100): # 100 timesteps
	xl = shift(x,1) # I'm going to compare each point on x to the point to its left and right
	xr = shift(x,-1) # so this function makes an array of each of those points
	for j,val in enumerate(xr): #nested loops :( # bad
		x[j] = (xl[j]^val) # it's faster if I use the val rather than index the array; can't figure out how to store both xl and xr in val
		print(j)
	if np.mod(i,9)==0:	
		plt.plot(np.array(list(map(int, x)))*2.5+i ,'x') # plot the output every 9 iterations
plt.ylabel('Iteration')

print("--- %s seconds ---" % (time.clock() - start_time)) # record how long it took. - fkn ages
plt.show()

#(if I comment out the plotting commands it's about 0.05s)
