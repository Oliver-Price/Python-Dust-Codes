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
        drawopts = "Data Region Enclosed"
        sav_bool = False
        test_mode = True
    elif (saveexists == True):
        with open(savefile, 'rb') as f:
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
                        "Data Region Enclosed","Dust Phase Angles","No Image"]
            drawopts_ans = easygui.choicebox(dmsg, dtitle, dchoices)
            if drawopts_ans != None:
                drawopts = drawopts_ans
        
        if reply == 'Toggle Data Save': sav_bool = not sav_bool       
        if reply == 'Toggle Test Mode': test_mode = not test_mode
        if reply == 'OK All': break  
            
    with open(savefile, 'wb') as f:
        pickle.dump([betau, betal, bno, simtu, simtl, tno, threshold,
                     drawopts, sav_bool, test_mode], f)
    
    return (betau, betal, bno, simtu, simtl, tno, drawopts, sav_bool, test_mode)

#%%

def simulation_setup_lorfrag(savefile):
    
    saveexists = os.path.exists(savefile)
    
    if (saveexists == False):
        betau = 1; betal = 0.1; bno = 50
        simtu = 10; simtl = 1; tno = 50
        bsta = 1; ftime = 50
        threshold = None #formerly used for phase space reduction, now defunct
        drawopts = "Data Region Enclosed"
        sav_bool = False
        test_mode = True
    elif (saveexists == True):
        with open(savefile, 'rb') as f:
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
            bsta = sparameters[10]
            ftime = sparameters[11]
            
    while 1:
        
        message = ('Beta From ' + str(betal) + ' to ' + str(betau)
        +' : ' + str(bno) + ' values\n')
        message += ('Ejection Times From ' + str(simtl) + ' to ' + str(simtu)
        +' : ' + str(tno) + ' values\n')
        message += ('Pre-Fragmentation Beta: ' + str(bsta) + ' \n')
        message += ('Exposure Time: ' + str(ftime) + ' \n')
        message += ('Draw Options: ' + drawopts + ' \n')
        message += ('Save Data: ' + str(sav_bool) + ' \n')
        message += ('Test Mode: ' + str(test_mode) + ' \n')
        reply = easygui.buttonbox(msg= message,
                                  title='Choosing Simulation Parameters',
              choices=('Change Beta Values', 'Change Ejection Time Values',
                       'Change Frag Beta', 'Change Frag Time',
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
        
        if reply == 'Change Frag Beta':
    
            fbmsg = 'Enter Pre-Fragmentation Beta'
            fbtitle = 'Change Pre-Fragmentation Beta'
            fbVal = easygui.enterbox(fbmsg,fbtitle)
            
            while 1:
                if fbVal == None: break #exit if cancel pressed
                errbfmsg = ""
                if fbVal.strip() == "": #check if entered
                    errbfmsg = ('Enter Pre-Fragmentation Beta')
                if errbfmsg == "": break
                fbVal = easygui.enterbox(errbfmsg,fbtitle)
                
            if fbVal != None:
                bsta = float(fbVal)
                
        if reply == 'Change Frag Time':
    
            ftmsg = 'Enter Pre-Fragmentation Beta'
            fttitle = 'Change Pre-Fragmentation Beta'
            ftVal = easygui.enterbox(ftmsg,fttitle)
            
            while 1:
                if ftVal == None: break #exit if cancel pressed
                errftmsg = ""
                if ftVal.strip() == "": #check if entered
                    errftmsg = ('Enter Pre-Fragmentation Beta')
                if errftmsg == "": break
                ftVal = easygui.enterbox(errftmsg,fttitle)
                
            if ftVal != None:
                ftime = float(ftVal)
                              
        if reply == 'Change Draw Options':

            dmsg = "Choose draw options"
            dtitle = "Choosing draw options"
            dchoices = ["Synchrones Only", "Syndynes only",
                        "Synchrones and Syndynes",
                        "Wide Spaced Synchrones and Syndynes",
                        "Synchrones, Syndynes and Data Points",
                        "Data Region Enclosed","Dust Phase Angles","No Image"]
            drawopts_ans = easygui.choicebox(dmsg, dtitle, dchoices)
            if drawopts_ans != None:
                drawopts = drawopts_ans
        
        if reply == 'Toggle Data Save': sav_bool = not sav_bool       
        if reply == 'Toggle Test Mode': test_mode = not test_mode
        if reply == 'OK All': break  
            
    with open(savefile, 'wb') as f:
        pickle.dump([betau, betal, bno, simtu, simtl, tno, threshold,
                     drawopts, sav_bool, test_mode, bsta, ftime], f)
    
    return (betau, betal, bno, simtu, simtl, tno, drawopts, sav_bool, test_mode, bsta, ftime)

#%%

def simulation_setup_FLM(savefile):
    
    saveexists = os.path.exists(savefile)
    stdict = {True: 'Test',False: 'Save'}
    
    if (saveexists == False):
        t = 6; bfrau = 1.5; bfral = 0.5; bno = 10
        h = 1.2; cu = 140; cl = 30; cno = 10
        s0u = 0.4; s0l = 0.2; s0no = 10; sc = 0.125
        test_mode = True
    elif (saveexists == True):
        with open(savefile, 'rb') as f:
            sparameters = pickle.load(f)
            t = sparameters[0]
            bfrau = sparameters[1]
            bfral = sparameters[2]
            bno = sparameters[3]
            h  = sparameters[4]
            cu = sparameters[5]
            cl = sparameters[6]
            cno = sparameters[7]
            s0u = sparameters[8]
            s0l = sparameters[9]
            s0no = sparameters[10]
            sc = sparameters[11]
            test_mode = sparameters[12]
            
    while 1:
        
        message = ('Ejection time: ' + str(t) + ' \n')
        message += ('Beta of fragments from ' + str(bfrau) + ' to ' + str(bfral)
        +' : ' + str(bno) + ' values\n')
        message += ('h value: ' + str(h) + ' \n')
        message += ('c from ' + str(cu) + ' to ' + str(cl)
        +' : ' + str(cno) + ' values\n')
        message += ('S0 from ' + str(s0u) + ' to ' + str(s0l)
        +' : ' + str(s0no) + ' values\n')
        message += ('Sc value: ' + str(sc) + ' \n')
        message += ('Test Mode: ' + stdict[test_mode] + ' \n')
        reply = easygui.buttonbox(msg= message,
                                  title='Choosing Simulation Parameters',
                       choices=('Set Ejection Time value', 'Set Fragment Beta Range',
                                'Set h Value', 'Set c Range',
                                'Set S0 Range', 'Set Sc Value',
                                'Toggle Save/Test Mode', 'OK All'),
                       image=None)
                       
        if reply == None: sys.exit()
        if reply == 'Set Fragment Beta Range':
    
            bmsg = 'Choose Fragment Beta Range'
            btitle = 'Set Fragment Beta Range'
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
                bfrau = float(bfieldValues[1])
                bfral = float(bfieldValues[0])
                bno = int(bfieldValues[2])
                
        if reply == 'Set c Range':
    
            tmsg = "Choose c Values"
            ttitle = "Changing c Values"
            tfieldNames = ["c Lower", "c Upper",
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
                cu = float(tfieldValues[1])
                cl = float(tfieldValues[0])
                cno = int(tfieldValues[2])
                
        if reply == 'Set S0 Range':
    
            smsg = "Choose S0 Values"
            stitle = "Changing S0 Values"
            sfieldNames = ["S0 Lower", "S0 Upper",
                          "Number of Datapoints"]
            sfieldValues = []  #values to be assigned
            sfieldValues = easygui.multenterbox(smsg,stitle, sfieldNames)
            
            while 1:
                if sfieldValues == None: break #exit if cancel pressed
                errsmsg = ""
                for i in range(len(sfieldNames)):
                    if sfieldValues[i].strip() == "": #check if entered
                                    errsmsg += ('"%s" is a required field.\n\n'
                                    % sfieldNames[i])
                if errsmsg == "": break
                sfieldValues = easygui.multenterbox(errsmsg, stitle,
                                                    sfieldNames, sfieldValues)
            if sfieldValues != None:
                s0u = float(sfieldValues[1])
                s0l = float(sfieldValues[0])
                s0no = int(sfieldValues[2]) 
                
        if reply == 'Set Ejection Time value':
    
            etmsg = 'Set Ejection Time value'
            ettitle = 'Set Ejection Time value'
            etVal = easygui.enterbox(etmsg,ettitle)
            
            while 1:
                if etVal == None: break #exit if cancel pressed
                erretmsg = ""
                if etVal.strip() == "": #check if entered
                    erretmsg = ('Set Ejection Time value')
                if erretmsg == "": break
                etVal = easygui.enterbox(erretmsg,ettitle)
                
            if etVal != None:
                t = float(etVal)

        if reply == 'Set h Value':
    
            hmsg = 'Set h Value'
            htitle = 'Set h Value'
            hVal = easygui.enterbox(hmsg,htitle)
            
            while 1:
                if hVal == None: break #exit if cancel pressed
                errhmsg = ""
                if hVal.strip() == "": #check if entered
                    errhmsg = ('Set h Value')
                if errhmsg == "": break
                hVal = easygui.enterbox(errhmsg,htitle)
                
            if hVal != None:
                h = float(hVal) 

        if reply == 'Set Sc Value':
    
            scmsg = 'Set Sc Value'
            sctitle = 'Set Sc Value'
            scVal = easygui.enterbox(scmsg,sctitle)
            
            while 1:
                if scVal == None: break #exit if cancel pressed
                errscmsg = ""
                if scVal.strip() == "": #check if entered
                    errscmsg = ('Set Sc Value')
                if errscmsg == "": break
                scVal = easygui.enterbox(errscmsg,sctitle)
                
            if scVal != None:
                sc = float(scVal) 
                
        if reply == 'Toggle Test/Save Mode': test_mode = not test_mode
        if reply == 'OK All': break  
            
    with open(savefile, 'wb') as f:
        pickle.dump([t,bfrau,bfral,bno,h,cu,cl,cno,s0u,s0l,s0no,sc,test_mode], f)
    
    return (t,bfrau,bfral,bno,h,cu,cl,cno,s0u,s0l,s0no,sc,test_mode)