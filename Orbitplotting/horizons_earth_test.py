# -*- coding: utf-8 -*-
import numpy as np
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

#%% FORWARD
csv_ecliptic = np.genfromtxt ('earth_hor.csv', delimiter=",")
csv_equator = np.genfromtxt ('earth_hor_frame.csv', delimiter=",")

theta = np.radians(23.4392911)

a = csv_equator[:,2]*np.cos(theta) + csv_equator[:,3]*np.sin(theta)
b = csv_equator[:,3]*np.cos(theta) - csv_equator[:,2]*np.sin(theta)