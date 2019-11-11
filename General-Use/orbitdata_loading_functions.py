#**************************
#ORBIT DATA LOADING METHODS
#**************************

import os
import numpy as np
import easygui
import requests
import datetime
import astropy.time
import pickle
import sys
from IPython import embed
from conversion_routines import mon2num

#%%

#***********
#VECTOR DATA
#***********

def orb_vector(denom, observer, savefolder, datafolder, horiz = '', opts = ''):
    
    #gives loading parameters for specific data type requested
    if "obs" in opts:
        title = 'Select astrometric observer orbit data file'
        msg = 'Format: <Observer>_orbit_<dates>_xyzvxvyvz'
        savename = 'data_' + denom + '_' + observer \
        + '_orbit_xyzvxvyvz'
    else:
        title = 'Select astrometric comet orbit data file'
        msg = 'Format: <Denom>_<dates>_xyzvxvyvz'
        savename = 'data_' + denom + '_cometorbit_xyzvxvyvz'        
    if "eq" in opts:
        savename += '_EQ'
        msg += '_EQ'
    if "lt" in opts:
        savename += '_LT'
        msg += '_LT'
    if "d10" in opts:
        savename += '_coarse'
        msg += '_coarse'
    savename += '.npy'
                
    #checks if save data exsits   
    savefile = os.path.join(savefolder, savename)
    if not os.path.exists(savefolder): os.makedirs(savefolder)
    saveexists = os.path.exists(savefile)
    
    if saveexists == False: #generate data if it doesn't
        
        loadmsg = "Choose where to load orbit data from:"
        loadchoices = ["Load from File","Import from Horizons"]
        reply = easygui.buttonbox(loadmsg, choices=loadchoices)

        if reply == "Import from Horizons":
            
            if horiz == '':
                sys.exit("Please supply JPL Horizons tag")
            
            #saving the end date ensures all orbit files match up
            datesavename = denom + '_enddate.pickle'
            datesavefile = os.path.join(savefolder, datesavename)
            datesaveexists = os.path.isfile(datesavefile)
            
            if datesaveexists == False:
                msg = "Enter End Date for data download"
                title = "Fetching Horizons data from batch"
                fieldNames = ["Day", "Month", "Year"]
                fieldLengths = [2,2,4] #constrains digits in numbers
                fieldValues = []  #values to be assigned
                fieldValues = easygui.multenterbox(msg,title, fieldNames)
                
                # conduct sanity check on dates
                while 1:
                    if fieldValues == None: break #exit if cancel pressed
                    errmsg = ""
                    for i in range(len(fieldNames)):
                        if fieldValues[i].strip() == "": #check if entered
                            errmsg += ('"%s" is a required field.\n\n'
                            % fieldNames[i])
                        if len(fieldValues[i].strip()) != fieldLengths[i]: #check length
                            errmsg +=  ('"%s" must be a %d digit number.\n\n'
                                % (fieldNames[i], fieldLengths[i]))

                    try:    #check date is a real date
                        newDate = datetime.datetime(int(fieldValues[2]),
                                                    int(fieldValues[1]),
                                                    int(fieldValues[0]))
                    except ValueError:
                        errmsg += ('Date must be real.\n\n')
                    if errmsg == "": break #if no problems found
                    fieldValues = easygui.multenterbox(errmsg, title,
                                                       fieldNames, fieldValues)
                print (fieldValues)
                
                with open(datesavefile, 'wb') as f:
                    pickle.dump(fieldValues, f) #save date
                    
            else:
                with open(datesavefile,'rb') as f:
                    fieldValues = pickle.load(f) #otherwise, load date
            
            print (fieldValues)
            #astropy date format makes it easier to find start date
            endtime =  astropy.time.Time(datetime.datetime(int(fieldValues[2]),
                       int(fieldValues[1]), int(fieldValues[0]), 0 , 0, 0))
            
            #url is for website request
            #textsavename is name of file on local disk            
            url = 'http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='                         
            
            if "obs" in opts: #choose between earth and comet
                if observer == 'Earth':
                    url += '\'399\'' #EARTH
                    textsavename = os.path.join(datafolder, 'Earth_orbit_')
                if observer == 'Venus':
                    url += '\'299\'' #EARTH
                    textsavename = os.path.join(datafolder, 'Venus_orbit_')
                if observer == 'Stereo_A':
                    url += '\'-234\'' #EARTH
                    textsavename = os.path.join(datafolder, 'Stereo_A_orbit_')
                if observer == 'Stereo_B':
                    url += '\'-235\'' #EARTH
                    textsavename = os.path.join(datafolder, 'Stereo_B_orbit_')
                if observer == 'Soho':
                    url += '\'-21\'' #EARTH
                    textsavename = os.path.join(datafolder, 'Soho_orbit_')
                if observer == 'ISS':
                    url += '\'-125544\'' #EARTH
                    textsavename = os.path.join(datafolder, 'ISS_orbit_')
                        
            else:
                if (horiz[0] != '\''): #this just ensures the formatting is all cleared up
                    url = url + '\''
                url += horiz.replace(' ', '%20') #COMET
                if (horiz[-1] != '\''):
                    url = url + '\''
                textsavename = os.path.join(datafolder, denom + '_') 
                
            url += '&MAKE_EPHEM=\'YES\'&TABLE_TYPE=\'VECTOR\'&CENTER=\'500@10\''
            
            if "d10" in opts: #choose between 1/10 minute data
                stepsiz = '\'10%20min\'' #10 minute dt for 600 days
                dta = astropy.time.TimeDelta(600, format = 'jd')
            else:
                stepsiz = '\'1%20min\'' #1 minute dt for 60 days
                dta = astropy.time.TimeDelta(60, format = 'jd') 
                
            starttime = endtime - dta
            url += '&START_TIME=\'' + starttime.isot[:10] + '\''
            url += '&STOP_TIME=\'' + endtime.isot[:10]  + '\''
            url += '&STEP_SIZE=' + stepsiz
            
            #clears up the format a little
            textsavename += (starttime.isot[:10].replace('-','_') + '_'
                            + endtime.isot[:10].replace('-','_')
                            + '_xyzvxvyvz')
                            
            if "eq" in opts:
                url += '&REF_PLANE=\'FRAME\''
                textsavename += '_EQ'
                
            url += '&OUT_UNITS=\'AU-D\'&VECT_TABLE=\'2\''
            
            if "lt" in opts:
                url += '&VECT_CORR=\'LT\''
                textsavename += '_LT'
            
            if "d10" in opts:
                textsavename += '_10'
            
            textsavename += '.txt'
            print(url)
            html = requests.get(url).text
            datastart = html.find('$$SOE') + 6 #selects relevant data
            dataend = html.find('$$EOE') - 1
            
            with open(textsavename, "w") as text_file:
                text_file.write(html[datastart:dataend]) #write to file         

            datapath = textsavename #tells it where to look
            
        if reply == "Load from File": #depreciated manual data loading mode
            #gets user to find file
            datalook = os.path.join(datafolder, '*.txt')   
            datapath = easygui.fileopenbox(msg=msg, title=title,
                                            default = datalook)
            
        with open(datapath, "r") as dat:
                
            #reads raw data as string 
            rawdata = dat.readlines()
            datasize = int(len(rawdata)/3)
            data = np.empty((datasize,13),dtype = float)
            
            print(rawdata)
            
            #converts string to numpy array
            for hrow in range(0,datasize,1):
                data[hrow,0] = rawdata[3*hrow][0:17]            #jd
                data[hrow,1] = rawdata[3*hrow][25:29]           #year
                data[hrow,2] = mon2num(rawdata[3*hrow][30:33])  #month
                data[hrow,3] = rawdata[3*hrow][34:36]           #day
                data[hrow,4] = rawdata[3*hrow][37:39]           #hour
                data[hrow,5] = rawdata[3*hrow][40:42]           #minute
                data[hrow,6] = rawdata[(hrow*3)+1][4:26]
                data[hrow,7] = rawdata[(hrow*3)+1][30:52]
                data[hrow,8] = rawdata[(hrow*3)+1][56:78]
                data[hrow,9] = rawdata[(hrow*3)+2][4:26]
                data[hrow,10] = rawdata[(hrow*3)+2][30:52]
                data[hrow,11] = rawdata[(hrow*3)+2][56:78]
                data[hrow,12] = np.linalg.norm(data[hrow,6:9])
            np.save(savefile,data)   #saves data for future loading

    elif saveexists == True: #if save exists, load that
        data = np.load(savefile)
        print ('Loading saved data')
    return data
    
#%%

#*************
#OBSERVER DATA
#*************

def orb_obs(denom, observer, savefolder, datafolder, horiz, idlmode = False, venusmode = False, mercurymode = False):

    #checks if save data exsits
    if venusmode == True:
        savename = 'data_' + denom + '_' + observer + '_venuscoords'
    elif mercurymode == True:
        savename = 'data_' + denom + '_' + observer + '_mercurycoords'    
    else:
        savename = 'data_' + denom + '_' + observer + '_celestialcoords'
    if idlmode == True:
        savename += '_idl'
    savename += '.npy'
    
    savefile = os.path.join(savefolder, savename)
    saveexists = os.path.isfile(savefile)

    if saveexists == False:
        
        loadmsg = "Choose where to load orbit data from:"
        loadchoices = ["Load from File","Import from Horizons"]
        reply = easygui.buttonbox(loadmsg, choices=loadchoices)

        if reply == "Import from Horizons":
            
            datesavename = denom + '_enddate.pickle'
            datesavefile = os.path.join(savefolder, datesavename)
            datesaveexists = os.path.isfile(datesavefile)

            if datesaveexists == False:
                msg = "Enter End Date for data download"
                title = "Fetching Horizons data from batch"
                fieldNames = ["Day", "Month", "Year"]
                fieldLengths = [2,2,4] #constrains digits in numbers
                fieldValues = []  #values to be assigned
                fieldValues = easygui.multenterbox(msg,title, fieldNames)
                
                # make sure fields are entered and correct length
                while 1:
                    if fieldValues == None: break #exit if cancel pressed
                    errmsg = ""
                    for i in range(len(fieldNames)):
                        if fieldValues[i].strip() == "": #check if entered
                            errmsg += ('"%s" is a required field.\n\n'
                            % fieldNames[i])
                        if len(fieldValues[i].strip()) != fieldLengths[i]: #check length
                            errmsg +=  ('"%s" must be a %d digit number.\n\n'
                                % (fieldNames[i], fieldLengths[i]))

                    try:
                        newDate = datetime.datetime(int(fieldValues[2]),
                                                    int(fieldValues[1]),
                                                    int(fieldValues[0]))
                    except ValueError:
                        errmsg += ('Date must be real.\n\n')
                    if errmsg == "": break #if no problems found
                    fieldValues = easygui.multenterbox(errmsg, title,
                                                       fieldNames, fieldValues)
                                                       
                with open(datesavefile, 'wb') as f:
                    pickle.dump(fieldValues, f)
                    
            else:
                with open(datesavefile, 'rb') as f:
                    fieldValues = pickle.load(f)
    
            endtime =  astropy.time.Time(datetime.datetime(int(fieldValues[2]),
                       int(fieldValues[1]), int(fieldValues[0]), 0 , 0, 0))
            
            url = 'http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='                         
            
            if venusmode == True:
                url += '\'299\''
                textsavename = os.path.join(datafolder, 'Venus_') 
            elif mercurymode == True:
                url += '\'199\''
                textsavename = os.path.join(datafolder, 'Mercury_')     
            
            else:
                if (horiz[0] != '\''): #this just ensures the formatting is all cleared up
                    url += '\''
                
                url +=  horiz.replace(' ', '%20')
                
                if (horiz[-1] != '\''):
                    url +=  '\''
           
                textsavename = os.path.join(datafolder, denom + '_') 
    
                
            url += '&MAKE_EPHEM=\'YES\'&TABLE_TYPE=\'OBS\'&CENTER='
            
            if observer == 'Earth':
                url += '\'500\'' #EARTH
                textsavename += 'from_Earth_'
            if observer == 'Stereo_A':
                url += '\'500@-234\'' #STEREO-A
                textsavename += 'from_Stereo_A_'
            if observer == 'Stereo_B':
                url += '\'500@-235\'' #STEREO-B
                textsavename += 'from_Stereo_B_'
            if observer == 'Soho':
                url += '\'500@-21\'' #EARTH
                textsavename += 'from_Soho_'
            if observer == 'ISS':
                    url += '\'500@-125544\'' #EARTH
                    textsavename = os.path.join(datafolder, 'ISS_orbit_')
                        
    
            starttime = endtime - astropy.time.TimeDelta(60, format = 'jd')
            url += '&START_TIME=\'' + starttime.isot[:10] + '\''
            url += '&STOP_TIME=\'' + endtime.isot[:10]  + '\''
            
            if (idlmode == False):
                url += '&STEP_SIZE=\'1%20min\'&QUANTITIES=\'1,28\''
                url +=  '&ANG_FORMAT=\'DEG\'&CSV_FORMAT=\'YES\''
                
                textsavename += (starttime.isot[:10].replace('-','_')
                             + '_' + endtime.isot[:10].replace('-','_')
                             + '_OBSERVER_PY.txt' )
            elif (idlmode == True):
                url += '&STEP_SIZE=\'1%20min\'&QUANTITIES=\'1,19,20,27\''
                url +=  '&CSV_FORMAT=\'YES\''
                
                textsavename += (starttime.isot[:10].replace('-','_')
                             + '_' + endtime.isot[:10].replace('-','_')
                             + '_OBSERVER_IDL.txt' )
                             
            print (url)
            print (textsavename)
            
            html = requests.get(url).text
            datastart = html.find('$$SOE') + 6 #selects relevant data
            dataend = html.find('$$EOE') - 1
            
            with open(textsavename, "w") as text_file:
                text_file.write(html[datastart:dataend])            

            datapath = textsavename    

        if reply == "Load from File": #gets user to find file
                
            title = 'Select observer orbit data file'
            msg = 'Format: <denom>_<dates>_OBSERVER_PY'
            datalook = os.path.join(datafolder, '*.txt')
            datapath = easygui.fileopenbox(msg=msg, title=title,default = datalook)

        if (idlmode == False):
            with open(datapath, "r") as dat:
                rawdata = dat.readlines()
                data = np.empty((len(rawdata),8),dtype = float)
                for row in range(0,len(rawdata),1):
                    data[row,0] = rawdata[row][1:5]             #year
                    data[row,1] = mon2num(rawdata[row][6:9])    #month
                    data[row,2] = rawdata[row][10:12]           #day
                    data[row,3] = rawdata[row][13:15]           #hour
                    data[row,4] = rawdata[row][16:18]           #minute
                    data[row,5] = rawdata[row][23:32]           #ra
                    data[row,6] = rawdata[row][33:42]           #dec
                    data[row,7] = rawdata[row][44:53]           #orbit plane angle
                np.save(savefile,data)    #saves data for future loading
                
        if (idlmode == True):
            
            data = 'void'
            
    elif saveexists == True: #if save exists, load that
        data = np.load(savefile)
        print ('Loading saved data')
    return data

#%%

#*************
#OBSERVER DATA
#*************

def orb_obs_extra(denom, observer, savefolder, datafolder):

    #checks if save data exsits
    savename = 'data_' + denom + '_' + observer + '_celestialapparents'
    savename += '.npy'
    
    savefile = os.path.join(savefolder, savename)
    saveexists = os.path.isfile(savefile)

    if saveexists == False:
         
        title = 'Select observer orbit data file'
        msg = 'Format: <denom>_<dates>_OBSERVER_PY_APP'
        datalook = os.path.join(datafolder, '*.txt')
        datapath = easygui.fileopenbox(msg=msg, title=title,default = datalook)

        with open(datapath, "r") as dat:
            rawdata = dat.readlines()
            data = np.empty((len(rawdata),11),dtype = float)
            for row in range(0,len(rawdata),1):
                data[row,0] = rawdata[row][1:5]             #year
                data[row,1] = mon2num(rawdata[row][6:9])    #month
                data[row,2] = rawdata[row][10:12]           #day
                data[row,3] = rawdata[row][13:15]           #hour
                data[row,4] = rawdata[row][16:18]           #minute
                data[row,5] = rawdata[row][23:32]           #ra
                data[row,6] = rawdata[row][33:42]           #dec
                data[row,7] = rawdata[row][44:53]           #ra_app
                data[row,8] = rawdata[row][54:63]           #dec_app
                data[row,9] = rawdata[row][64:73]
                data[row,10] = rawdata[row][74:82]
            np.save(savefile,data)    #saves data for future loading
            
    elif saveexists == True: #if save exists, load that
        data = np.load(savefile)
        print ('Loading saved data')
    return data

#%%

#*************
#NEW METHODS
#*************

def orb_vector_new(denom, observer, savefolder, datafolder, opts = ''):
    
    #gives loading parameters for specific data type requested
    if "obs" in opts:
        savename = 'data_' + denom + '_' + observer \
        + '_orbit_xyzvxvyvz'
    else:
        savename = 'data_' + denom + '_cometorbit_xyzvxvyvz'        
    if "eq" in opts:
        savename += '_EQ'
    if "d10" in opts:
        savename += '_coarse'
    savename += '.npy'
                
    #checks if save data exsits   
    savefile = os.path.join(savefolder, savename)
    if not os.path.exists(savefolder): os.makedirs(savefolder)
    saveexists = os.path.exists(savefile)
    
    if saveexists == False: #generate data if it doesn't
        sys.exit("Please load data via Jupyter notebook")
        
    elif saveexists == True: #if save exists, load that
        data = np.load(savefile)
        print ('Loading saved data')
    return data

def orb_obs_new(denom, observer, savefolder, datafolder, opts = ''):

    #checks if save data exsits
    if 'v' in opts:
        savename = 'data_' + denom + '_' + observer + '_venuscoords'
    if 'm' in opts:
        savename = 'data_' + denom + '_' + observer + '_mercurycoords'    
    else:
        savename = 'data_' + denom + '_' + observer + '_celestialcoords'
    savename += '.npy'
    
    savefile = os.path.join(savefolder, savename)
    saveexists = os.path.isfile(savefile)

    if saveexists == False:
        sys.exit("Please load data via Jupyter notebook")
        
    elif saveexists == True: #if save exists, load that
        data = np.load(savefile)
        print ('Loading saved data')
        
    return data