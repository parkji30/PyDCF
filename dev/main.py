import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from dcf_python import *

# plt.rcParams.update({'font.size': 16})

## Find location ##
os.chdir("/Users/a16472/Desktop/MDCF/dev")


## Loading the polarization angle data
data = fits.open("L1M10_0.1.fits")[0].data


## Coordinates and Size ##
y_cen = (450)
x_cen = (100)
rad = 50

# Taking a smaller region from the entire map.
data_region = data_cut(x_cen, y_cen, rad, data, show=False)

# Calculating the structure function analysis of the smaller region.
dr, dphi = cos_disp_calculations(data_region, ds_scale=1)


## Corrected Turbulent to Ordered Ratio
fit0_01 = 20
fitf_01 = 25
uncorrected_turbulent_ratio, turb_cof = MDCF_fit(dr, dphi, 'L1M10_0.2', edge_length=1.0, \
                                         beam_res=0.2, fit0=fit0_01, fitf=fitf_01, beam=True)

N = turbulent_cells(turb_cof, 1.45, 0.2/2.35)
print("Turbulent Cells:", round(N))

print("Corrected Turbulent to Ordered Ratio")
print("------------------------------------")
if N > 1:
    print(np.sqrt(round(N) * uncorrected_turbulent_ratio))
else:
    print(np.sqrt(uncorrected_turbulent_ratio))