#default parameters
import easygui
import pickle
import os
import sys

def simulation_setup(savefile):
    
    saveexists = os.path.exists(savefile)
    
    if (saveexists == False):
        betau = 1; betal = 0.1; bno = 50
        simtu = 10; simtl = 1; tno = 50
        threshold = None #formerly used for phase space reduction, now defunct
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
            threshold = sparameters[6]
            drawopts = sparameters[7]
            sav_bool = sparameters[8]
            test_mode = sparameters[9]
            
    while 1:
        
        message = ('Beta From ' + str(betal) + ' to ' + str(betau)
        +' : ' + str(bno) + ' values\n')
        message += ('Ejection Times From ' + str(simtl) + ' to ' + str(simtu)
        +' : ' + str(tno) + ' values\n')
        message += ('Draw Options: ' + drawopts + ' \n')
        message += ('Save Data: ' + str(sav_bool) + ' \n')
        message += ('Test Mode: ' + str(test_mode) + ' \n')
        reply = easygui.buttonbox(msg= message,
                                  title='Choosing Simulation Parameters',
              choices=('Change Beta Values', 'Change Ejection Time Values',
                       'Change Draw Options', 'Toggle Data Save',
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
                                       
        if reply == 'Change Draw Options':

            dmsg = "Choose draw options"
            dtitle = "Choosing draw options"
            dchoices = ["Synchrones Only", "Syndynes only",
                        "Synchrones and Syndynes",
                        "Wide Spaced Synchrones and Syndynes",
                        "Synchrones, Syndynes and Data Points",
                        "Data Region Enclosed","No Image"]
            drawopts_ans = easygui.choicebox(dmsg, dtitle, dchoices)
            if drawopts_ans != None:
                drawopts = drawopts_ans
        
        if reply == 'Toggle Data Save': sav_bool = not sav_bool       
        if reply == 'Toggle Test Mode': test_mode = not test_mode
        if reply == 'OK All': break  
            
    with open(savefile, 'w') as f:
        pickle.dump([betau, betal, bno, simtu, simtl, tno, threshold,
                     drawopts, sav_bool, test_mode], f)
    
    return (betau, betal, bno, simtu, simtl, tno, drawopts, sav_bool, test_mode)
    