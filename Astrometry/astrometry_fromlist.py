import requests
import urllib
import json 
import time
import sys
import os

png_location = r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\pngs"


session_login = requests.post('http://nova.astrometry.net/api/login',
                  data={'request-json': json.dumps({"apikey": "rqleqqsuqhatwjni"})})
                  
session_key = session_login.text.split("\"")[-2]
                 
image_list = os.listdir(png_location)
    
for image_no in xrange(1, len(image_list)):

    print image_no    
    
    image_png_url = ("http://www.mssl.ucl.ac.uk/~op2/" + image_list[image_no])
    
    image_submit = requests.post("http://nova.astrometry.net/api/url_upload",
                  data={'request-json': json.dumps({"session": session_key,
                  "url": image_png_url})})
                 
    image_sub_id = image_submit.text.split(":")[2][1:8]   
    image_sub_stat_url = r"http://nova.astrometry.net/api/submissions/" + image_sub_id
    image_sub_stat = requests.get(image_sub_stat_url).text
    
    test = True
    while test == True:
        image_sub_stat = requests.get(image_sub_stat_url).text
        time.sleep(0.5)
        test = ("null" in image_sub_stat) or ("None" in image_sub_stat)
    
    image_saveas = os.path.join(r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B_recalibrated",
                                image_list[image_no] + ".fits")   
    image_fits_file_loc = (r"http://nova.astrometry.net/new_fits_file/" +
                            image_sub_stat.split("jobs")[1][4:11])

    urllib.urlretrieve(image_fits_file_loc, image_saveas)
    while os.path.getsize(image_saveas) < 1000000:
        urllib.urlretrieve(image_fits_file_loc, image_saveas)
        time.sleep(0.5)
        
        
#{"session": "nlybgrbj2p3bgi3ankmqy9593i5sftb1",
# "url": "http://imageupper.com/s03/1/2/W14663519672416591_1.png"}
#urllib.urlretrieve(test, "test.fits")
#{"session": "o3eqg241tgcxuvllvxfikzcuek44ahtm", "url": "http://imageupper.com/s02/1/3/A14663520302419729_42.png"}