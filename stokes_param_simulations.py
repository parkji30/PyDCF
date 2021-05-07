import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy
import os
import aplpy as ap

os.chdir("C:/Users/16472/Desktop/Queens 2020 Fall/msc. Thesis/MAGNETIC_FIELD_DISPERSION/images")


## Functions
def psi(Q, U):
    """
    Calculate azimuth angle.
    """

    return np.arctan2(Q, U)/2


def Imshow(image, **kwargs):
    """
    Simple function to an image.
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(image)
    plt.colorbar()
    plt.show()

def coterminal_angle

## Calculation Pol Angle

q = fits.open("VelaC_500_intermediate_regrid_30as_pix_Q.fits")[0].data
u = fits.open("VelaC_500_intermediate_regrid_30as_pix_U.fits")[0].data
masking = fits.open("VelaC_deepmap_mask.fits")[0].data
true_masking = np.load("true_masking.npy")

psi_ang = psi(q, u)[:220][:]

psi_ang = psi_ang *true_masking
Imshow(psi_ang, origin='lower')

## Lauras

psi_map = fits.open("VelaC_500_intermediate_regrid_30as_pix_ang.fits")
angular_map = psi_map[0].data
plt.figure(figsize=(16, 12))
plt.title("Laura's Angular Map")
plt.imshow(angular_map)
plt.colorbar()
plt.show()


## Now Compare the two

Imshow(angular_map - psi_ang)


## For simulation
q = fits.open("ccldL10B5M10n0025v2fwhm5_polQ.fits")[0].data
u = fits.open("ccldL10B5M10n0025v2fwhm5_polU.fits")[0].data

psi_ang = psi(q, u)
Imshow(psi_ang)

# Shift it by 45 deg, or 0.78539

Imshow(psi_ang)