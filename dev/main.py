import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit

from dcf_python import *
plt.rcParams.update({'font.size': 16})

os.chdir("/Users/a16472/Desktop/dcf_python/testing")

data = fits.open("L10_0.65.fits")[0].data
data = data[::6, ::6]
y_cen = 260 / 6
x_cen = 180 / 6
rad = 25

dt = data_cut(x_cen, y_cen, rad, data, show=True)

dr, dphi = cos_disp_calculations(dt, ds_scale=1)
multi_fit(dr, dphi, 'Region', ds_scale=1, outer_distance=21, fit0=13, fitf=21, show=True)