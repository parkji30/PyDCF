import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit

from dcf_python import *
plt.rcParams.update({'font.size': 16})

## Find location ##
os.chdir("/Users/a16472/Desktop/dcf_python/testing")

## Coordinates and Size ##
y_cen = (260) / 2
x_cen = (160) / 2

# y_cen = (512*3/4+25) /2
# x_cen = (512/4-25) /2

# y_cen = 100
# x_cen = 325
rad = 100


data = fits.open("L10_0.2.fits")[0].data
data = data[::2, ::2]

dt = data_cut(x_cen, y_cen, rad, data, show=True)
dr, dphi = cos_disp_calculations(dt, ds_scale=2)

fit0_01 = 20
fitf_01 = 25
uncorrected_turbulent_ratio, turb_cof  = single_fit(dr, dphi, 'L1M10_0.2', edge_length=1.0, \
                                         beam_res=0.2, fit0=fit0_01, fitf=fitf_01, beam=True)

N = turbulent_cells(turb_cof, 1.45, 0.2/2.35)
print("Turbulent Cells:", round(N))

print("Corrected Turbulent to Ordered Ratio")
print("------------------------------------")
if N > 1:
    print(np.sqrt(round(N) * uncorrected_turbulent_ratio))
else:
    print(np.sqrt(uncorrected_turbulent_ratio))