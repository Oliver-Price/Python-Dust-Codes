#####################################
# UPLOADS A BATCH OF IMAGES TO ASTROMETRY.NET FROM FTP SERVER
#####################################

import requests
import json 
import time
import os
import numpy as np
import urllib
import pickle

orig_location = r"C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\notext"
save_loc = r"C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\fits"

session_login = requests.post('http://nova.astrometry.net/api/login',
                  data={'request-json': json.dumps({"apikey": "rqleqqsuqhatwjni"})})
                  
session_key = session_login.text.split("\"")[-2]

orig_list = os.listdir(orig_location)

todo_list = orig_list[293:]

#%%
'''
badsave = r"C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\badlist2.pickle"
donesave = r"C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\donelist2.pickle"

with open (badsave, 'rb') as fp:
    bad_list = pickle.load(fp)
    
with open (donesave, 'rb') as fp:
    done_list = pickle.load(fp)
    
not_list = bad_list + done_list
todo_list = np.setdiff1d(np.array(orig_list),np.array(not_list)).tolist()
'''
#%%

out_list = []
url_list = []
submit_list = []
image_sub_ids = []
image_sub_stat_urls = []
image_sub_stats = []
image_fits_file_locs = []

#%%

start = 0; end = len(todo_list)
for image_no in range(start,end):

    out_list.append(os.path.join(save_loc,('.').join(todo_list[image_no].split('.')[:-1])+'.fits'))

    #fits online
    url_list.append('http://www.mssl.ucl.ac.uk/~op2/Cometpics/' + todo_list[image_no])

    print (todo_list[image_no])   
        
    submit_list.append(requests.post("http://nova.astrometry.net/api/url_upload",
                      data={'request-json': json.dumps({"session": session_key,
                      "url": url_list[image_no-start], "downsample_factor": 1})}))
                     
    image_sub_ids.append(submit_list[image_no-start].text.split(":")[2][1:8])  
    image_sub_stat_urls.append(r"http://nova.astrometry.net/api/submissions/" + image_sub_ids[image_no-start])
    image_sub_stats.append(requests.get(image_sub_stat_urls[image_no-start]).text)
    
    print(submit_list[image_no-start].text)
    print('\n')
    time.sleep(5)
    
#%%
    
start = 0; end = len(todo_list)

for image_no in range(start,end):#len(todo_list)):
    print('retrieving: ' + todo_list[image_no])
    
    if "}" not in image_sub_stat_urls[image_no-start]:
    
        while 1:
            if (len(requests.get(image_sub_stat_urls[image_no-start]).json()['jobs']) == 0):
                break
            
            job_no = requests.get(image_sub_stat_urls[image_no-start]).json()['jobs'][0]
            image_sub_stat = requests.get('http://nova.astrometry.net/api/jobs/' + str(job_no)).text
            time.sleep(2)
            if ("success" in image_sub_stat) or ("failure" in image_sub_stat):
                break
        
        if ("success" in image_sub_stat):
            image_fits_file_locs = r"http://nova.astrometry.net/new_fits_file/" + str(job_no)
            print("fetching fits")
            urllib.request.urlretrieve(image_fits_file_locs,out_list[image_no-start])
            done_list.append(orig_list[image_no])
        else:
            print("Image was no good")
            bad_list.append(orig_list[image_no])
            
    else:
        print("Image was no good")
        bad_list.append(orig_list[image_no])
        
    print('\n')

#%%

with open(badsave, 'wb') as fp:
    pickle.dump(bad_list, fp)
    
with open(donesave, 'wb') as fp:
    pickle.dump(done_list, fp)

#%%

#{"session": "nlybgrbj2p3bgi3ankmqy9593i5sftb1",
# "url": "http://imageupper.com/s03/1/2/W14663519672416591_1.png"}
#urllib.urlretrieve(test, "test.fits")
#{"session": "o3eqg241tgcxuvllvxfikzcuek44ahtm", "url": "http://imageupper.com/s02/1/3/A14663520302419729_42.png"}