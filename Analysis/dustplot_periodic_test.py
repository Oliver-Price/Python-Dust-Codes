#dustplot FITS analysis tool
import easygui
import sys
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import string

#choosing comet data to use
inputfilefolder = "C:\PhD\Comet_data\Input_files\*pt1.txt"
inputfile = easygui.fileopenbox(default = inputfilefolder)

#reading main comet parameters
with open(inputfile, "r") as c:
    cdata = c.readlines()
    comname = cdata[30][12:]
    comdenom = cdata[31][13:-2]
    imagedir = cdata[24][18:-2]
    orbitdir = cdata[25][23:-2]
    idlsav = cdata[26][25:-2]
    pysav = cdata[27][24:-2]
    obsloc = cdata[34][19:24]
    horiztag = cdata[40][10:]
    
dustplotsave = os.path.join(imagedir, 'dustplots')
if not os.path.exists(dustplotsave):
    sys.exit("No Dustplots Calculated")
    
dustplot_in = easygui.fileopenbox(default =
                                    os.path.join(dustplotsave, r'*.fits'))
hdulist = fits.open(dustplot_in)
dustarr = (hdulist[0].data)
dustarr = dustarr.astype(float)

simresdir = os.path.join(pysav, 'simres')
fitsinfile = os.path.basename(dustplot_in)
filebase = fitsinfile[:string.find(fitsinfile,'filter')][:-2]
simin = os.path.join(simresdir, filebase + '.npy')
simres = np.load(simin)

zerolocs = np.where(dustarr == 0)
for n in xrange(0,np.shape(zerolocs)[1]):
    dustarr[zerolocs[0][n],zerolocs[1][n]] = np.nan
    
dust_shape_1 = np.shape(dustarr)[1]
dust_x_sect_1 = np.empty((dust_shape_1,2))

for d0 in xrange(0, dust_shape_1):
    dust_x_sect_1[d0,0] = d0 
    dust_x_sect_1[d0,1] = np.nanmean(dustarr[:,d0])

av_vals = dust_x_sect_1[:,1][~np.isnan(dust_x_sect_1[:,1])]
t_vals = simres[:,0,0][~np.isnan(dust_x_sect_1[:,1])]
fft_vals = np.fft.fft(av_vals).real

a = np.sin(np.arange(45)*1.0/18* np.pi)
b = np.arange(45)
c = np.fft.fft(a)

matplotlib.rcParams.update({'font.size': 35})
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
tes = ax.plot(t_vals, av_vals)
#ts = ax.plot(np.arange(np.size(t_vals)),fft_vals)
plt.xlabel('Time since Ejection (days)')
plt.ylabel('Average filtered pixel value of synchrone')




    