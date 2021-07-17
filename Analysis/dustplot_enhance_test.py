# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

de = srcolors[:,:,1]-gaussian_filter(srcolors[:,:,1], sigma=20)

#%%
num_bins = 5
n, bins, patches = plt.hist(de, num_bins, facecolor='blue', alpha=0.5)

#%%
plt.imshow(np.rot90(np.clip(de,-5,5).T,2))