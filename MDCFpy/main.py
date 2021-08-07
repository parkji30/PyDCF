import numpy as np
import os
os.chdir("/Users/a16472/Desktop/MDCF/MDCFpy")
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from dcf_python import *
from fitting_tools import *

## Loading the polarization angle data
data = fits.open("L1M10_0.1.fits")[0].data


## Coordinates and Size ##
y_cen = (450)
x_cen = (100)
rad = 50

cloud_depth = 1.56

# Taking a smaller region from the entire map.
data_region = data_cut(x_cen, y_cen, rad, data, show=False)


fit0_01 = 20
fitf_01 = 25

# Calculating the structure function analysis of the smaller region.
uncorrected_turbulent_ratio, turbulent_coefficient = MDCF_fit(data_region,
                                                              pixel_scale = 10/512,
                                                              edge_length = 0.5,
                                                              beam_size = 0.1,
                                                              fit0 = 20,
                                                              fitf = 25)



N = turbulent_cells(turbulent_coefficient, cloud_depth, 0.1/2.35)
print("Turbulent Cells:", round(N))

print("\nCorrected Turbulent to Ordered Ratio")
print("------------------------------------")
if N > 1:
    print(np.sqrt(round(N) * uncorrected_turbulent_ratio))
else:
    print(np.sqrt(uncorrected_turbulent_ratio))
