import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
os.chdir("/Users/a16472/desktop/mdcf/PyDCF")

from PyDCF import *
from dispersion_analysis import *

## Loading the polarization angle data
data = fits.open("L1M10_0.1.fits")[0].data
velocity = fits.open("L1M10_sigmav_0.1.fits")[0].data
print(np.mean(velocity))
density = fits.open("L1M10_meanrho_0.1.fits")[0].data


## Region Snipping
y_cen = (280)
x_cen = (140)
rad = 80

# Taking a smaller region from the entire map.
data_pol_region = data_cut(x_cen, y_cen, rad, data, show=True)
data_v_region = data_cut(x_cen, y_cen, rad, velocity, show=False)
data_rho_region = data_cut(x_cen, y_cen, rad, density, show=False)


## Declare PyDCF
'''
beam resolution and pixel scale must be in the same units.
'''

pold1 = PyDCF(data_pol_region,
              data_v_region,
              data_rho_region,
              beam_resolution = 0.1,
              pixel_scale = 10/512)


pold1.calculate_angular_dispersions()

## Fit
pold1.HH09_fit(18, 25, 1.51)


## Calculate Magnetic Field Strength
print(str(pold1.ClassicalDCF_calculation()*1e6) + " Microgauss")
print(str(pold1.SkalidisDCF_calculation()*1e6) + " Microgauss")
print(str(pold1.HH09DCF_calculation()*1e6) + " Microgauss")
