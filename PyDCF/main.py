import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
os.chdir("/Users/a16472/desktop/mdcf/PyDCF")

from PyDCF import *
from dispersion_analysis import *

## Loading the polarization angle data
data = fits.open("L1M10_0.1.fits")[0].data
velocity = fits.open("L1M10_meanrho_0.1.fits")[0].data
density = fits.open("L1M10_meanrho_0.1.fits")[0].data


## Region Snipping
y_cen = (260)
x_cen = (140)
rad = 40

# Taking a smaller region from the entire map.
data_pol_region = data_cut(x_cen, y_cen, rad, data, show=False)
data_v_region = data_cut(x_cen, y_cen, rad, velocity, show=False)
data_rho_region = data_cut(x_cen, y_cen, rad, density, show=False)


## Declare PyDCF


pold1 = PyDCF(data_pol_region, data_v_region, data_rho_region, 0.1, 10/512)
pold1.calculate_angular_dispersions()
pold1.HH09_fit(18, 25, 1.51)


## Calculate Magnetic Field Strength
pold1.HH09DCF_calculation()
pold1.classicalDCF_calculation()