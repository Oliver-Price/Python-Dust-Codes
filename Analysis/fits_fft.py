# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from astropy.io import fits

betau = 2.5
betal = 0.1
bno = 400
simtu = 6.0
simtl = 1.0
tno = 400

datafolder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Linux Data\SOHO Blue\fits'
image_list = sorted(os.listdir(datafolder))
image_total = len(image_list)

image_id = 36

image_procname = image_list[image_id].split('.')[0][:-7]
image_mapped = os.path.join(datafolder,image_list[image_id])
print(image_procname)

hdulist = fits.open(image_mapped)
srcolors = (hdulist[0].data).astype(float)
srcolors = np.swapaxes(srcolors,0,1)

tvals = np.linspace(simtl, simtu, tno)
bvals = np.logspace(np.log10(betal), np.log10(betau), num=(bno))
#%%

betaval = 1.0

row = abs(bvals-betaval).argmin()
rowS = np.where(srcolors[:,row]<=5000)[0][0]
rowE = np.where(srcolors[:,row]>=400)[0][-1]


#%%
a = tvals[rowS:rowE]
b = srcolors[rowS:rowE,row]

# Number of samplepoints
N = rowE-rowS
# sample spacing
T = 1.0 / N
x = np.linspace(0.0, N*T, N)
bf = scipy.fftpack.fft(b)
xf = np.nan_to_num(np.reciprocal(np.linspace(0.0, N/2, N/2)))[1:]
xf = xf*N*(tvals[1]-tvals[0])*24
yf = 2.0/N * np.abs(bf[:N//2])[1:]
  
fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1)
plottitle = image_procname[:15] + "  beta = " + str(betaval)

plt.tight_layout(h_pad=2.0,rect=[0, 0.03, 1, 0.95])
plt.suptitle(plottitle, fontsize=16)
ax1.semilogy(tvals[rowS:rowE],srcolors[rowS:rowE,row])
ax1.set_xlabel('Age (Days)')
ax1.set_ylabel('Brightness')

ax2.semilogx(xf,yf)
ax2.set_xlabel('Period of features(Hours)')
ax2.set_ylabel('Relative Power')

plt.show()
fig.subplots_adjust(bottom=0.15,left=0.15)

#plotsavename = os.path.join(plotsave,type + "_" + image_procname[:15] + "_beta_" + str(betaval).replace('.','\'') + '.png')
#fig.savefig(plotsavename,dpi=100)
