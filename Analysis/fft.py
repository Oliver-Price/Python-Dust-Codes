# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

betau = 3.5
betal = 0.2
bno = 800
simtu = 20.0
simtl = 2.0
tno = 800

type = 'MGN'

datafolder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Linux Data\HI-1'
simresfolder = os.path.join(datafolder,'simres')
srcolorsfolder = os.path.join(datafolder,'srcolors-' + type)
plotsave = os.path.join(datafolder,'plots')

image_list = sorted(os.listdir(srcolorsfolder))
image_total = len(image_list)

image_id = 25

image_procname = image_list[image_id].split('.')[0][:-7]
image_mapped = os.path.join(srcolorsfolder,image_procname + "_mapped.npy")
image_simressave = os.path.join(simresfolder, image_procname+ ".npy")

print(image_procname)

simres = np.load(image_simressave)
srcolors = np.load(image_mapped)

tvals = np.linspace(simtl, simtu, tno)
bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))
#%%

betaval = 1.0

row = abs(bvals-betaval).argmin()
rowS = np.where(srcolors[:,row,0]==1)[0][0]
rowE = np.where(srcolors[:,row,0]==1)[0][-1]


#%%
a = tvals[rowS:rowE]
b = srcolors[rowS:rowE,row,1]

# Number of samplepoints
N = rowE-rowS
# sample spacing
T = 1.0 / N
x = np.linspace(0.0, N*T, N)
bf = scipy.fftpack.fft(b)
xf = np.nan_to_num(np.reciprocal(np.linspace(0.0, N/2, N/2)))[1:]
xf = xf*N*(tvals[1]-tvals[0])*24
yf = 2.0/N * np.abs(bf[:N//2])[1:]
  
fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1,sharex = False)
plottitle = type + " " + image_procname[:15] + "  beta = " + str(betaval)

plt.tight_layout(h_pad=2.0,rect=[0, 0.03, 1, 0.95])
plt.suptitle(plottitle, fontsize=16)
ax1.plot(tvals[rowS:rowE],srcolors[rowS:rowE,row,1])
ax1.set_xlabel('Age (Days)')
ax1.set_ylabel('Brightness')

ax2.semilogx(xf,yf)
ax2.set_xlabel('Period of features(Hours)')
ax2.set_ylabel('Relative Power')

plt.show()
fig.subplots_adjust(bottom=0.15,left=0.15)

plotsavename = os.path.join(plotsave,type + "_" + image_procname[:15] + "_beta_" + str(betaval).replace('.','\'') + '.png')
fig.savefig(plotsavename,dpi=100)
