#EASY PREP FOR LOADING ORBIT DATA FOR IMAGETIMEHEADERS
sys.path.append(r"C:\PhD\Python\Python-Dust-Codes\General-Use")  
from orbitdata_loading_functions import orb_vector, orb_obs
import easygui
import numpy as np
import os
import sys

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
    obslocstr = cdata[34][19:]
    horiztag = cdata[40][10:]

#choose observer locations
bool_locs = np.array([(('EARTH' in obslocstr) or ('Earth' in obslocstr)),
                 (('STEREO-A' in obslocstr) or ('Stereo-A' in obslocstr)
                 or ('STEREO_A' in obslocstr) or ('Stereo_A' in obslocstr)),
                 (('STEREO-B' in obslocstr) or ('Stereo-B' in obslocstr)
                 or ('STEREO_B' in obslocstr) or ('Stereo_B' in obslocstr)),
                 (('SOHO' in obslocstr) or ('Soho' in obslocstr))])
name_locs = np.array(['Earth', 'Stereo_A', 'Stereo_B', 'Soho'])
case_locs = np.size(np.nonzero(bool_locs))
if case_locs > 1:
    obsmsg = "Please select observer location"
    obschoices = name_locs[bool_locs].tolist()
    obsloc = easygui.buttonbox(obsmsg, choices=obschoices)
    imagedir = os.path.join(imagedir, obsloc)
elif case_locs == 1:
    obsloc = name_locs[bool_locs][0]
else: sys.exit("No Good Observer Location")

obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,lt')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'lt')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag, idlmode = True)
