# -*- coding: utf-8 -*-

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.cluster import DBSCAN

file_base = r"C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\Small" #r before the " is to ignore string formatting commands e.g \n

im_out_dir = os.path.join(file_base,'dbscan')
if not os.path.exists(im_out_dir):
    os.makedirs(im_out_dir)

img_list = os.listdir(file_base)
imgs = [s for s in img_list if "." in s]

eps = 150
msm = 1050

overlay = True

subdir = "e" + str(eps) + "m" + str(msm) + "o" + str(overlay)

im_out_dir = os.path.join(im_out_dir,subdir)
if not os.path.exists(im_out_dir):
    os.makedirs(im_out_dir)

for i in range(len(imgs)):
    
    file_in = os.path.join(file_base,imgs[i])
    image = sp.misc.imread(file_in)
    
    #remove alpha channel and flatten
    image_reduced = image[..., :3]
    img_flat = sp.reshape(image_reduced, (image_reduced.shape[0]*image_reduced.shape[1],3)).astype(float)

    dbf = DBSCAN(eps = eps, min_samples = msm)
    dbimage_indices_flat = dbf.fit_predict(img_flat) #get the clustered data
    dbimage_indices = sp.reshape(dbimage_indices_flat, (image.shape[0], image.shape[1]))
    
    if overlay ==  True:  
        
        db_overlay = np.copy(image_reduced)
        
        #find the background with mode and color it red
        db_com = sp.stats.mode(dbimage_indices,axis=None)[0][0]
        
        if db_com == dbimage_indices.min():
            db_overlay[np.where(dbimage_indices == dbimage_indices.max())] = np.array([255,0,0])
    
        if db_com == dbimage_indices.max():
            db_overlay[np.where(dbimage_indices == dbimage_indices.min())] = np.array([255,0,0])
    
        dbsave = db_overlay
        
    else:
        
        dbsave = dbimage_indices       
    
    #blow up and save
    db_large_save = sp.misc.imresize(dbsave,(dbsave.shape[0]*5,dbsave.shape[1]*5,3))
    file_save = os.path.join(im_out_dir,imgs[i])
    sp.misc.imsave(file_save,db_large_save)