import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from regions import PixCoord, RectanglePixelRegion

def Imshow(image, **kwargs):
    """
    Simple function display to an image.
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(image, origin='lower')
    plt.colorbar()
    plt.show()


def plt_show(x, y):
    """
    Show relation between x and y.
    """
    plt.figure(figsize = (10, 10))
    plt.plot(x, y)
    plt.show()


def hist_show(array, title='', xlabel='Value'):
    """
    Histogram of array.
    """
    plt.figure()
    plt.hist(array, bins=350)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Counts")
    plt.show()


def report_stats(array):
    """
    Basic Stats Calculatiiosn
    """
    array = array[~np.isnan(array)]
    mean = np.mean(array)
    std = np.std(array)
    median = np.median(array)
    print("Mean: ", mean)
    print("Median: ", median)
    print("Standard Deviation: ", std)


def total_gauss_function(x, a, W, delta):
    """
    Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2*(delta**2 + 2*W**2)))


def gauss_function(x, a, sigma):
    """
    Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2 * sigma**2))


def beam_gauss_function(x, a):
    """
    Beam Gaussian fit function.
    """
    W = 0.56626
    return a * np.exp(-(x)**2 / (2*(2 * W**2)))


def linear_fit(x, m, b):
    """
    Linear Fit Function.
    """
    return m*x + b


def data_cut(x_cen, y_cen, rad, image, show=False):
    if show:
        fig, ax = plt.subplots(figsize=(10,6))
        region = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
        plt.imshow(image, cmap='hsv', vmin=-np.pi/2, vmax=np.pi/2)
        plt.colorbar()
        region.plot(ax=ax, color='white')
        plt.show()
    reg = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
    mask = reg.to_mask()
    mask = reg.to_mask(mode='center')
    dt = mask.cutout(image)
    return dt