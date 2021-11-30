import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits

def linear_fit(x, m, b):
    """
    Linear Fit Function.
    """
    return m * x + b


def gauss_function(x, a, sigma):
    """
    Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2 * sigma**2))


def total_gauss_function(x, a, W, delta):
    """
    Turbulent Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2*(delta**2 + 2*W**2)))


def large_scale_fit(l, b_ratio_sq, delta, a):
    """
    """
    W = 0.05
    cloud_depth = 1.46
    return np.sqrt(np.pi) * (b_ratio_sq) * ( delta**3 / ((delta**2 + 2*W**2) * cloud_depth)) \
            * ( 1 - np.exp(-l**2/(2*(delta**2 + 2*W**2)))) + (a * l**2)


def turbulent_autocorrelation(distance, b_ratio, delta, W, ED):
    """
    b_strength: turbulent to large scale magnetic field ratio
    delta: delta

    Used to fit the gaussian for the 3rd plot.
    """
    return  np.sqrt(2*np.pi) * b_ratio**2 * (delta**3 / ((delta**2 + 2*W**2) * ED)) \
            * (np.exp(-distance**2 / (2*(delta**2 + 2*W**2))))


def turbulent_cells(delta, cloud_dep, beam_res):
    """
    """
    return ((delta**2 + 2*(beam_res**2)) * cloud_dep) / (np.sqrt(2*np.pi) * delta**3)


def skalidis_DCF(field_density, sigma_v, sigma_pol):
    return np.sqrt(2*np.pi*field_density) * sigma_v / np.sqrt(sigma_pol)


def classical_DCF(field_density, sigma_v, sigma_pol):
    return np.sqrt(4*np.pi*field_density) * sigma_v / sigma_pol


def modified_DCF(field_density, velocity_los, b_ratio):
    """
    Return magnetic field of desired region.
    """
    return np.sqrt(4*np.pi*field_density) * velocity_los * b_ratio**(-0.5)