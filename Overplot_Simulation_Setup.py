#default parameters
import easygui
import pickle
import os
import sys

def simulation_setup(savefile):
    
    saveexists = os.path.exists(savefile)
    
    if (saveexists == False):
        betau = 1; betal = 0.001; bno = 50
        simtu = 100; simtl = 1; tno = 50; tspace = 'Linear'
        threshold = 10
        drawopts = "Synchrones Only"
        sav_bool = False
        test_mode = True
    elif (saveexists == True):
        with open(savefile) as f:
            sparameters = pickle.load(f)
            betau = sparameters[0]
            betal = sparameters[1]
            bno = sparameters[2]
            simtu = sparameters[3]
            simtl = sparameters[4]
            tno = sparameters[5]
            tspace = sparameters[6]
            threshold = sparameters[7]
            drawopts = sparameters[8]
            sav_bool = sparameters[9]
            test_mode = sparameters[10]
            
    while 1:
        
        message = ('Beta From ' + str(betal) + ' to ' + str(betau)
        +' : ' + str(bno) + ' values\n')
        message += ('Ejection Times From ' + str(simtl) + ' to ' + str(simtu)
        +' : ' + str(tno) + ' values\n')
        message += ('Ejection Time Spacing: ' + tspace + ' \n')
        message += ('Draw Options: ' + drawopts + ' \n')
        message += ('Blocking Threshold = ' + str(threshold)+ ' \n')
        message += ('Save Data: ' + str(sav_bool) + ' \n')
        message += ('Test Mode: ' + str(test_mode) + ' \n')
        reply = easygui.buttonbox(msg= message,
                                  title='Choosing Simulation Parameters',
              choices=('Change Beta Values', 'Change Ejection Time Values',
                       'Toggle Time Spacing', 'Change Draw Options',
                       'Change Threshold', 'Toggle Data Save',
                       'Toggle Test Mode', 'OK All'),
                       image=None)
                       
        if reply == None: sys.exit()
        if reply == 'Change Beta Values':
    
            bmsg = "Choose Beta Values"
            btitle = "Changing Beta Values"
            bfieldNames = ["Beta Lower", "Beta Upper",
                          "Number of Datapoints"]
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
                betau = float(bfieldValues[1])
                betal = float(bfieldValues[0])
                bno = int(bfieldValues[2])
                
        if reply == 'Change Ejection Time Values':
    
            tmsg = "Choose Ejection Time Values"
            ttitle = "Changing Ejection Time Values"
            tfieldNames = ["Ejection Time Lower", "Ejection Time Upper",
                          "Number of Datapoints"]
            tfieldValues = []  #values to be assigned
            tfieldValues = easygui.multenterbox(tmsg,ttitle, tfieldNames)
            
            while 1:
                if tfieldValues == None: break #exit if cancel pressed
                errtmsg = ""
                for i in range(len(tfieldNames)):
                    if tfieldValues[i].strip() == "": #check if entered
                                    errtmsg += ('"%s" is a required field.\n\n'
                                    % tfieldNames[i])
                if errtmsg == "": break
                tfieldValues = easygui.multenterbox(errtmsg, ttitle,
                                                    tfieldNames, tfieldValues)
            if tfieldValues != None:
                simtu = float(tfieldValues[1])
                simtl = float(tfieldValues[0])
                tno = int(tfieldValues[2])
                                       
        if reply == 'Toggle Time Spacing':
           
            if tspace == 'Logarithmic': tspace = 'Linear'                  
            elif tspace == 'Linear': tspace = 'Logarithmic' 
                                      
        if reply == 'Change Draw Options':

            dmsg = "Choose draw options"
            dtitle = "Choosing draw options"
            dchoices = ["Synchrones Only", "Syndynes only",
                        "Synchrones and Syndynes",
                        "Synchrones, Syndynes and Data Points",
                        "Data Region Enclosed","No Image"]
            drawopts_ans = easygui.choicebox(dmsg, dtitle, dchoices)
            if drawopts_ans != None:
                drawopts = drawopts_ans
                
        if reply == 'Change Threshold':   
          
            hmsg = "Choose Threshold for ignoring near-coma data"
            htitle = "Choosing cut-off Threshold"
            hreply = easygui.enterbox(hmsg,htitle)
            while 1:
                if hreply == None: break #exit if cancel pressed                 
                errmsg = ""
                try:
                    threshold = int(hreply)
                except ValueError:
                        errmsg += ("Threshold must be an Integer")
                if errmsg == "": break
                hreply = easygui.enterbox(errmsg, htitle)
            threshold = int(hreply)
        
        if reply == 'Toggle Data Save': sav_bool = not sav_bool       
        if reply == 'Toggle Test Mode': test_mode = not test_mode
        if reply == 'OK All': break  
            
    with open(savefile, 'w') as f:
        pickle.dump([betau, betal, bno, simtu, simtl, tno, tspace, threshold,
                     drawopts, sav_bool, test_mode], f)
    
    return (betau, betal, bno, simtu, simtl, tno, tspace, threshold, drawopts,
            sav_bool, test_mode)
            
    