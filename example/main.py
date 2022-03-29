import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import PyDCF

## Loading the polarization angle data
def data_cut(x_cen, y_cen, rad, image, show=False):
    '''
    This function is used to cut a square region from a map based
    on the radius and coordinates provided.

    The last argument, show, can be toggled to display the region.

    @type x_cen: Int
    @type y_cen: Int
    @type rad: Float
    @type image: Numpy Array
    @type show: Boolean
    @rtype: Numpy Array (Snipped Region)
    '''
    if show:
        fig, ax = plt.subplots(figsize=(10,6))
        region = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
        plt.imshow(image, cmap='hsv')
        plt.colorbar()
        region.plot(ax=ax, color='white')
        plt.show()
    reg = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
    mask = reg.to_mask()
    mask = reg.to_mask(mode='center')
    dt = mask.cutout(image)
    return dt

data = fits.open("L1M10_0.1.fits")[0].data
velocity = fits.open("L1M10_sigmav_0.1.fits")[0].data
density = fits.open("L1M10_meanrho_0.1.fits")[0].data

## Region Snipping
y_cen = (280)
x_cen = (140)
rad = 60

# Taking a smaller region from the entire map.
data_pol_region = data_cut(x_cen, y_cen, rad, data, show=True)
data_v_region = data_cut(x_cen, y_cen, rad, velocity, show=False)
data_rho_region = data_cut(x_cen, y_cen, rad, density, show=False)


## Declare PyDCF
'''
beam resolution and pixel scale must be in the same units.
'''
pold1 = PyDCF(polarization = data_pol_region,
              velocity = data_v_region,
              density = data_rho_region,
              beam_resolution = 0.1,
              pixel_scale = 10/512)


pold1.calculate_angular_dispersions()


## MDCF Fit
pold1.HH09_fit(fit0 = 18, fitf = 25, cloud_depth = 1.51)
pold1.HH09_parameters()


## Calculate Magnetic Field Strength
print(str(pold1.ClassicalDCF_calculation()*1e6) + " Microgauss")
print(str(pold1.SkalidisDCF_calculation()*1e6) + " Microgauss")
print(str(pold1.HH09DCF_calculation()*1e6) + " Microgauss")
x = pold1.correction_factors(10/1e6)