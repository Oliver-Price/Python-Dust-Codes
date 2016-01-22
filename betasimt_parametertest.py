#code for testing beta - t relationships

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import easygui
savefilefolder = r"C:\PhD\Python\Save_Data\c"
simin = easygui.fileopenbox(default = savefilefolder)
with open(simin) as f:
        sparameters = pickle.load(f)
        simres = sparameters[0]
        tmax = sparameters[1]
        bmax = sparameters[2]
        tno = sparameters[3]
        bno = sparameters[4]
        
tsiz = np.shape(simres)[0]
bsiz = np.shape(simres)[1]

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
        data[a,0] = simres[t,b,0] #simt
        data[a,1] = simres[t,b,1] #beta
        data[a,2] = abs(simres[t,b,12] - cimgral)/pwid #final RA in image widths
        data[a,3] = abs(simres[t,b,13] - cimgdecl)/phei #final DEC in image widths
        w = np.ceil(np.floor(data[a,2])/1000)
        h = np.ceil(np.floor(data[a,3])/1000)
        data[a,4] = toowide[w,h]
        data[a,5] = abs(simres[t,b,12] - simres[t-1,b,12])
        data[a,6] = abs(simres[t,b,13] - simres[t-1,b,13])
        data[a,7] = abs(simres[t,b,12] - simres[t,b-1,12])
        data[a,8] = abs(simres[t,b,13] - simres[t,b-1,13])
        e = np.ceil(np.floor(data[a,5])/1000000)
        f = np.ceil(np.floor(data[a,6])/1000000)
        g = np.ceil(np.floor(data[a,7])/1000000)
        h = np.ceil(np.floor(data[a,8])/1000000)
        data[a,9] = smaldel[e,f,g,h]
        data[a,10] = abs(simres[t,b,12] - comral)
        data[a,11] = abs(simres[t,b,13] - comdel)
        i = np.ceil(np.floor(data[a,10]/10)/1000000)
        j = np.ceil(np.floor(data[a,11]/10)/1000000)
        data[a,12] = smaldelcom[i,j]
        data[a,13] = -data[a,4] + data[a,9] + data[a,12]
        a +=1

#%%

simtl = simres[0,0,0]; simtu = simres[tno-1,0,0]
betal = simres[0,0,1]; betau = simres[0,bno-1,1]

decades = np.logspace(-4,4,9)
tdecl = decades[np.searchsorted(decades,simtl, side = 'right')-1]
tdecu = decades[np.searchsorted(decades,simtu)]
bdecl = decades[np.searchsorted(decades,betal, side = 'right')-1]
bdecu = decades[np.searchsorted(decades,betau)]

matplotlib.rcParams.update({'font.size': 35})
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(data[:,0], data[:,1],c=data[:,13],s=50, lw = 0)
plt.axis([tdecl, tdecu, bdecl, bdecu])
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel('Beta')
plt.xlabel('Time since emission (Days)')

#line, = plt.plot([12,60,80],[1,0.024,0.01], 'b')
