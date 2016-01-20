#code for testing beta - t relationships

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
savefile = r'C:\PhD\Python\Save_Data\simres'
simreso = np.load(savefile)
tsiz = np.shape(simreso)[0]
bsiz = np.shape(simreso)[1]

cimgral = 961.67533907229995
cimgdecl = 341.79721063243437
phei = 800
pwid = 1065

data = np.empty(((tsiz-1)*(bsiz-1),14),dtype = float)
a = 0
toowide = np.array([[0, 1], [1, 1]])
smaldel = np.ones((2,2,2,2))
smaldel[1,1,:,:] = 0
smaldel[:,:,1,1] = 0
comra = 3.503584069564437
comral = 957.14194152708217
comde = 0.2305931178802338
comdel = 335.58580423057583
smaldelcom = np.array([[1, 0], [0, 0]])

for t in xrange(1,tsiz):
    for b in xrange(1,bsiz):
        data[a,0] = simreso[t,b,0] #simt
        data[a,1] = simreso[t,b,1] #beta
        data[a,2] = abs(simreso[t,b,12] - cimgral)/pwid #final RA in image widths
        data[a,3] = abs(simreso[t,b,13] - cimgdecl)/phei #final DEC in image widths
        w = np.ceil(np.floor(data[a,2])/1000)
        h = np.ceil(np.floor(data[a,3])/1000)
        data[a,4] = toowide[w,h]
        data[a,5] = abs(simreso[t,b,12] - simreso[t-1,b,12])
        data[a,6] = abs(simreso[t,b,13] - simreso[t-1,b,13])
        data[a,7] = abs(simreso[t,b,12] - simreso[t,b-1,12])
        data[a,8] = abs(simreso[t,b,13] - simreso[t,b-1,13])
        e = np.ceil(np.floor(data[a,5])/1000000)
        f = np.ceil(np.floor(data[a,6])/1000000)
        g = np.ceil(np.floor(data[a,7])/1000000)
        h = np.ceil(np.floor(data[a,8])/1000000)
        data[a,9] = smaldel[e,f,g,h]
        data[a,10] = abs(simreso[t,b,12] - comral)
        data[a,11] = abs(simreso[t,b,13] - comdel)
        i = np.ceil(np.floor(data[a,10]/10)/1000000)
        j = np.ceil(np.floor(data[a,11]/10)/1000000)
        data[a,12] = smaldelcom[i,j]
        data[a,13] = -data[a,4] + data[a,9] + data[a,12]
        a +=1

#%%
matplotlib.rcParams.update({'font.size': 35})
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(data[:,0], data[:,1],c=data[:,13],s=300)
#plt.axis([0.001, 100, 0.01, 1000])
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel('Beta')
plt.xlabel('Time since emission (Days)')

#line, = plt.plot([12,60,80],[1,0.024,0.01], 'b')
