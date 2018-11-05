# -*- coding: utf-8 -*-
import numpy as np

p = np.linspace(-0.2, 0.2, 1000)
xv,yv = np.meshgrid(p,p)

xy = np.zeros((xv.size,3))
xy[:,0] = xv.reshape(xv.size)
xy[:,1] = yv.reshape(yv.size)