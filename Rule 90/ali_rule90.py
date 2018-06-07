import numpy as np
import matplotlib.pyplot as plt
import time

start=time.clock()

width = 100
time_steps = 100

curr_line = np.zeros(width,dtype=bool)
curr_line[np.uint(width/4)]=1
curr_line[np.uint(3*width/4)]=1
next_line = np.zeros(width,dtype=bool)
history = np.zeros(shape=(time_steps,width),dtype=bool)

history[0,:]=curr_line

for step in range(time_steps-1):
    for i in range(1,width-1):
        next_line[i] = history[step,i-1] != history[step,i+1]
    history[step+1,:] = next_line

end = time.clock()

print('Time elapsed: ', end-start)

fig,ax = plt.subplots(1)

ax.imshow(history,cmap = 'Greys')
plt.show(block=True)
