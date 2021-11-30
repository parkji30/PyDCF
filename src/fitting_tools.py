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


def data_cut(x_cen, y_cen, rad, image, show=False):
    '''
    Data cut version 1 which
    '''
    if show:
        fig, ax = plt.subplots(figsize=(10,6))
        region = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
        plt.imshow(image, cmap='hsv')#, vmin=-np.pi/2, vmax=np.pi/2)
        plt.colorbar()
        region.plot(ax=ax, color='white')
        plt.show()
    reg = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
    mask = reg.to_mask()
    mask = reg.to_mask(mode='center')
    dt = mask.cutout(image)
    return dt


def data_cut2(x_cen, y_cen, rad, image, show=False, vmin=0, vmax=10000, title='Default'):
    '''
    '''
    if show:
        fig, ax = plt.subplots(figsize=(10,6))
        region = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
        plt.title(title)
        plt.imshow(image, vmin=vmin, vmax=vmax)#, vmin=-np.pi/2, vmax=np.pi/2)
        plt.colorbar()
        region.plot(ax=ax, color='white')
        plt.show()
    reg = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
    mask = reg.to_mask()
    mask = reg.to_mask(mode='center')
    dt = mask.cutout(image)
    return dt
