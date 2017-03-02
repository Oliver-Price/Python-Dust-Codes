import requests
import urllib
import json 
import time
import sys
import os
import numpy as np

orig_location = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\Others\background_subtracted"
save_loc = r"C:\PhD\Comet_data\Comet_Lovejoy_C2011W3\Gallery\LASCO C3 Clear\Exptime_Sorted\Others\background_calibrated"

session_login = requests.post('http://nova.astrometry.net/api/login',
                  data={'request-json': json.dumps({"apikey": "rqleqqsuqhatwjni"})})
                  
session_key = session_login.text.split("\"")[-2]
                 
orig_list = os.listdir(orig_location)
done_list = os.listdir(save_loc)
todo_list = np.setdiff1d(np.array(orig_list),np.array(done_list)).tolist()
out_list = []
url_list = []
submit_list = []
image_sub_ids = []
image_sub_stat_urls = []
image_sub_stats = []
image_fits_file_locs = []

start = 0
for image_no in xrange(start,len(todo_list)):

    out_list.append(os.path.join(save_loc,orig_list[image_no]))

    #fits online
    url_list.append('ftp://ftp.mssl.ucl.ac.uk/pub/planet/op2/' + orig_list[image_no])

    #png online
    #url_list.append('ftp://ftp.mssl.ucl.ac.uk/pub/planet/op2/' + orig_list[0][:-5] + '_25_200.png')
    print orig_list[image_no-start]     
        
    submit_list.append(requests.post("http://nova.astrometry.net/api/url_upload",
                      data={'request-json': json.dumps({"session": session_key,
                      "url": url_list[image_no-start], "downsample_factor": 1, "center_ra": 263, "center_dec": -23, "radius": 15})}))
                     
    image_sub_ids.append(submit_list[image_no-start].text.split(":")[2][1:8])  
    image_sub_stat_urls.append(r"http://nova.astrometry.net/api/submissions/" + image_sub_ids[image_no-start])
    image_sub_stats.append(requests.get(image_sub_stat_urls[image_no-start]).text)
    
    print submit_list[image_no-start].text
    print '\n'
    time.sleep(0.5)
        
for image_no in xrange(start,len(todo_list)):
    print 'retrieving: ' + orig_list[image_no-start] 
    while 1:
        job_no = requests.get(image_sub_stat_urls[image_no-start]).json()['jobs'][0]
        image_sub_stat = requests.get('http://nova.astrometry.net/api/jobs/' + str(job_no)).text
        time.sleep(0.2)
        if ("success" in image_sub_stat) or ("failure" in image_sub_stat):
            break

    if ("success" in image_sub_stat):
        image_fits_file_locs = r"http://nova.astrometry.net/new_fits_file/" + str(job_no)
        urllib.urlretrieve(image_fits_file_locs, out_list[image_no-start])
    
    if ("failure" in image_sub_stat):
        print "Image was no good"
        
    print '\n'

       
#{"session": "nlybgrbj2p3bgi3ankmqy9593i5sftb1",
# "url": "http://imageupper.com/s03/1/2/W14663519672416591_1.png"}
#urllib.urlretrieve(test, "test.fits")
#{"session": "o3eqg241tgcxuvllvxfikzcuek44ahtm", "url": "http://imageupper.com/s02/1/3/A14663520302419729_42.png"}