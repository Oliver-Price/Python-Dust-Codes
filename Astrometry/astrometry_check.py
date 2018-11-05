# -*- coding: utf-8 -*-
################################
#BATCH CHECK ASTROMETRY.NET UPLOADS BY DOWNLOADING ANNOTATED DISPLAY
################################
import requests
import shutil
import os
imgout = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_A\HI-1\redo-astrometry\rg_check'

for image_no in range(1771559,1771584): # 1759850
    out_url = r'http://nova.astrometry.net/user_images/' + str(image_no) + '#redgreen'
    p = requests.get(out_url).text
    c = '/annotated_display/'
    if p.find(c) and p.find('fts') > 0:
        i_url = r'http://nova.astrometry.net/red_green_image_display/' + str(p[p.find(c)+19:p.find(c)+26])
        response = requests.get(i_url, stream=True)
        iout = os.path.join(imgout,p[p.find('fts')-22:p.find('fts')]+'jpg')
        with open(iout, 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
        del response