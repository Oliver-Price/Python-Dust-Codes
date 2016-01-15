#EASY PREP FOR LOADING ORBIT DATA FOR IMAGETIMEHEADERS

from orbitdata_loading_functions import orb_vector, orb_obs
import easygui

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
    
obsveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'obs,lt')
comveceq = orb_vector(comdenom, obsloc, pysav, orbitdir,
                      horiztag, opts = 'lt')
comobs = orb_obs(comdenom, obsloc, pysav, orbitdir, horiztag, idlmode = True)
