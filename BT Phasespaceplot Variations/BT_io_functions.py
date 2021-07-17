#default parameters
import easygui
import pickle
import os
import sys
import astropy
import datetime

#isot[0:16].replace('T',' at ')

def fixed_image_times(savefile, simtime_first, simtime_last):
    
    saveexists = os.path.exists(savefile)
    
    if (saveexists == False):
        img_start_time = astropy.time.Time(2457754.5,format='jd')
        img_end_time = astropy.time.Time(2457755.5,format='jd')
        beta_l = 0.1
        beta_u = 3
    elif (saveexists == True):
        with open(savefile, 'rb') as f:
            sparameters = pickle.load(f)
            img_start_time = sparameters[0]
            img_end_time = sparameters[1]   
            beta_l = sparameters[2]
            beta_u = sparameters[3]   
            
    while 1:
        
        message = ('Please choose start and end times\n\n')
        message += ('Observation Range:\n' +
                    simtime_first.isot[0:16].replace('T',' at ') + '\n' +
                    simtime_last.isot[0:16].replace('T',' at ') + '\n\n')
        message += ('Observation Range:\n' +
                    img_start_time.isot[0:16].replace('T',' at ') + '\n' +
                    img_end_time.isot[0:16].replace('T',' at ')+ '\n\n')
        message += ('Beta Range:\n' +
                   'Beta Lower: ' + str(beta_l) + '\n' +
                    'Beta Upper: ' + str(beta_u))
        
        reply = easygui.buttonbox(msg= message, title='Choose Image Timeframe',
                                  choices=('Change Start Time','Change End Time','Change Beta','OK All'),
                                  image=None)
                       
        if reply == None: sys.exit()
        if reply == 'Change Start Time':
            img_start_time = dategetter()
        if reply == 'Change End Time':
            img_end_time = dategetter()
        if reply == 'Change Beta':
            
            bmsg = "Choose Beta Values"
            btitle = "Changing Beta Values"
            bfieldNames = ["Beta Lower", "Beta Upper"]
            bfieldValues = []  #values to be assigned
            bfieldValues = easygui.multenterbox(bmsg,btitle, bfieldNames)
            
            while 1:
                if bfieldValues == None: break #exit if cancel pressed
                errbmsg = ""
                for i in range(len(bfieldNames)):
                    if bfieldValues[i].strip() == "": #check if entered
                                    errbmsg += ('"%s" is a required field.\n\n'
                                    % bfieldNames[i])
                if errbmsg == "": break
                bfieldValues = easygui.multenterbox(errbmsg, btitle,
                                                    bfieldNames, bfieldValues)
            
            if bfieldValues != None:                                       
                beta_u = float(bfieldValues[1])
                beta_l = float(bfieldValues[0])
            
        if reply == 'OK All': break  
            
    with open(savefile, 'wb') as f:
        pickle.dump([img_start_time,img_end_time,beta_l,beta_u], f)
    
    return (img_start_time,img_end_time,beta_l,beta_u)

#%%dategetter for things
def dategetter():
    
    img_t_msg = "Enter Date and Time of Image"
    img_t_title = "User input of image time"
    img_t_fieldNames = ["Year", "Month", "Day"]
    img_t_fieldLengths = [4,2,2] #constrains digits in numbers
    img_t_fieldValues = []  #values to be assigned
    img_t_fieldValues = easygui.multenterbox(img_t_msg, img_t_title,
                                             img_t_fieldNames)
    
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
                                        int(img_t_fieldValues[2]))
        except ValueError:
            errmsg += ('Date and Time must be real.\n\n')
        if errmsg == "": break #if no problems found
        img_t_fieldValues = easygui.multenterbox(errmsg,img_t_title,
                                           img_t_fieldNames,
                                           img_t_fieldValues)
    cday = int(img_t_fieldValues[2])
    cmonth = int(img_t_fieldValues[1])
    cyear = int(img_t_fieldValues[0])
    ctime = astropy.time.Time(datetime.datetime(cyear, cmonth, cday,
                                                0,0,0)) 
                                              
    return ctime