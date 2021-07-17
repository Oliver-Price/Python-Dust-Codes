#****************************
#Methods to get time of image
#****************************
import os
from scipy.io.idl import readsav
import datetime
import astropy.time
import easygui

#%% locate imgtimehdr save data from IDL

def image_time_yudish(comdenom,fitsinfile,idlsav):

            idlsavpath = os.path.join(idlsav,'Imagetimeheaders_savefile')
            idlsavpath = os.path.join(idlsavpath, comdenom)
            idlsavpath = os.path.join(idlsavpath, fitsinfile[len(comdenom)+1:-5])
            idlsavpath = idlsavpath + '_timeinfo.sav'
            idls = readsav(idlsavpath)
            
            #get time of comet in image and make an astropy time reference
            chour = int(idls.optocentre_time_str[0][0:2])
            cmin = int(idls.optocentre_time_str[0][3:5])
            cday = int(idls.optocentre_date[0][0:2])
            cmonth = int(idls.optocentre_date[0][3:5])
            cyear = int(idls.optocentre_date[0][6:11])
            ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                        chour , cmin, 0))
                                                        
            uncertainty_range_exists = True
            
            return idls,ctime,uncertainty_range_exists
            
#%% User inputs image time
            
def image_time_user():
    
    img_t_msg = "Enter Date and Time of Image"
    img_t_title = "User input of image time"
    img_t_fieldNames = ["Year", "Month", "Day", "Hour", "Minute"]
    img_t_fieldLengths = [4,2,2,2,2] #constrains digits in numbers
    img_t_fieldValues = []  #values to be assigned
    img_t_fieldValues = easygui.multenterbox(img_t_msg, img_t_title,
                                             img_t_fieldNames)
    # conduct sanity check on dates
    while 1:
        if img_t_fieldValues == None: break #exit if cancel pressed
        errmsg = ""
        for i in range(len(img_t_fieldNames)):
            if img_t_fieldValues[i].strip() == "": #check if entered
                errmsg += ('"%s" is a required field.\n\n'
                % img_t_fieldNames[i])
            if len(img_t_fieldValues[i].strip()) > img_t_fieldLengths[i]: #check length
                errmsg +=  ('"%s" must be at most a %d digit number.\n\n'
                    % (img_t_fieldNames[i], img_t_fieldLengths[i]))
    
        try:    #check date is a real date
            newDate = datetime.datetime(int(img_t_fieldValues[0]),
                                        int(img_t_fieldValues[1]),
                                        int(img_t_fieldValues[2]),
                                        int(img_t_fieldValues[3]),
                                        int(img_t_fieldValues[4]))
        except ValueError:
            errmsg += ('Date and Time must be real.\n\n')
        if errmsg == "": break #if no problems found
        fieldValues = easygui.multenterbox(errmsg,img_t_title,
                                           img_t_fieldNames,
                                           img_t_fieldValues)
                                           
    cmin = int(img_t_fieldValues[4])
    chour = int(img_t_fieldValues[3])
    cday = int(img_t_fieldValues[2])
    cmonth = int(img_t_fieldValues[1])
    cyear = int(img_t_fieldValues[0])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                chour , cmin, 0))                                                  
    uncertainty_range_exists = False
    
    return ctime,uncertainty_range_exists
    
#%% Image time from stereo filename
    
def image_time_stereo(filebase):
    
    csec = int(filebase[13:15])
    cmin = int(filebase[11:13])
    chour = int(filebase[9:11])
    cday = int(filebase[6:8])
    cmonth = int(filebase[4:6])
    cyear = int(filebase[0:4])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                chour , cmin, csec))
    uncertainty_range_exists = False                                          
    return ctime,uncertainty_range_exists

#%%Image time from yudish/earth
def image_time_filename_yuds(filebase):
    
    csec = 0
    cmin = int(filebase[21:23])
    chour = int(filebase[19:21])
    cday = int(filebase[16:18])
    cmonth = int(filebase[13:15])
    cyear = int(filebase[8:12])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                chour , cmin, csec))
    
    uncertainty_range_exists = False
    return ctime,uncertainty_range_exists

#%%Image time from yudish/earth
def image_time_filename_denom_compact(filebase):
    
    spl = filebase.split('_')
    csec = 0
    cmin = int(spl[2][2:4])
    chour = int(spl[2][0:2])
    cday = int(spl[1][6:])
    cmonth = int(spl[1][4:6])
    cyear = int(spl[1][0:4])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                chour , cmin, csec))
    
    uncertainty_range_exists = False
    return ctime,uncertainty_range_exists
