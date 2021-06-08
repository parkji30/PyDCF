import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits
from regions import PixCoord, CirclePixelRegion, RectanglePixelRegion
from scipy.optimize import curve_fit
from fitting_tools import *
plt.rcParams.update({'font.size': 16})

def calc_rel_angle_crossn(angle1, angle2, no_rescale=False):

    angle1 = np.array(angle1)
    angle2 = np.array(angle2)

    n = len(angle1)

    if n == 1:

        x1 = (-1.0) * np.sin(angle1[0])
        y1 = np.cos(angle1[0])
        x2 = (-1.0) * np.sin(angle2[0])
        y2 = np.cos(angle2[0])
        v1 = np.array([x1, y1, 0])
        v2 = np.array([x2, y2, 0])
        C = np.cross(v1, v2)

        CdC = np.dot(C, C)
        vdgr = np.dot(v1, v2)
        d_ang0 = np.arctan2(np.sqrt(CdC), vdgr)

        return np.array([d_ang0])

    elif n > 1:

        x1 = (-1.0) * np.sin(angle1.reshape(1, n))
        y1 = np.cos(angle1.reshape(1, n))
        x2 = (-1.0) * np.sin(angle2.reshape(1, n))
        y2 = np.cos(angle2.reshape(1, n))
        v1 = np.array([x1, y1, np.zeros((1, n))])
        v2 = np.array([x2, y2, np.zeros((1, n))])
        vi = np.asmatrix(v1).T
        vf = np.asmatrix(v2).T

        try:
            C = np.cross(vi, vf)
        except:
            print("crossing error!")

        CdC = np.sum(C * C, 1)

        vdgr = []
        for i in range(len(vi)):
            vector = v1[0][0][i] * v2[0][0][i] + \
                        v1[1][0][i] * v2[1][0][i] + \
                        v1[2][0][i] * v2[2][0][i]
            vdgr.append(vector)

        vdgr = np.array(vdgr)
        d_ang0 = np.arctan2(np.sqrt(CdC), np.abs(vdgr))
        return d_ang0


def cos_disp_calculations(data, ds_scale):
    """
    """
    x, y, pix_ang, dphi = [], [], [], []
    c = 0

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not np.isnan(data[i][j]):
                x.append(i)
                y.append(j)
                pix_ang.append(data[i][j])

    x = np.array(x)
    y = np.array(y)
    ang = np.array(pix_ang)

    nump = len(ang) 
    W = 2.5 / 2.35 # arc seconds
    delta_r = []
    delta_phi = []
    phi = ang

    for i in range(nump):
        delta_x_arr = x[i] - x[(i+1):(nump)]
        delta_y_arr = y[i] - y[(i+1):(nump)]
        delta_r_arr = np.sqrt(delta_x_arr**2 + delta_y_arr**2)
        sz_phi = len(delta_x_arr)
        phi_ref = np.repeat(phi[i], sz_phi)

        if len(phi_ref) > 0:
            delta_phi_arr = calc_rel_angle_crossn(phi_ref, phi[(i+1):(nump)])

        delta_r.append(delta_r_arr)
        delta_phi.append(delta_phi_arr)

    delta_r = np.array(delta_r)
    delta_phi = np.array(delta_phi[:-1]) # Last value is added twice for some reason.

    delta_r = np.concatenate(delta_r).ravel() * 10 / 512 * ds_scale # CONVERT THIS TO UNITS OF PARSEC
    delta_phi = np.concatenate(delta_phi).ravel()
    return delta_r, delta_phi


def multi_fit(delta_r, delta_phi, ttl, ds_scale, outer_distance, fit0=7, fitf=17, show=False):
    """
    """
    bin_range = (np.linspace(0, outer_distance, 21) + 0.5) * 10 / 512 

    # Binned Statistics calculation for the turbulent to ordered ratio.
    cos_disp, bin_edges_cos, bin_number_cos = stats.binned_statistic(delta_r, np.cos(delta_phi), 'mean', bins = bin_range)
    cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic(delta_r**2, np.cos(delta_phi), 'mean', bins = bin_range**2)
    cos_disp, bin_edges, bin_number_cos = stats.binned_statistic(delta_r, np.cos(delta_phi), 'mean', bins = bin_range)

    bin_edges_sq = np.insert(bin_edges_sq, 0, 0)
    cos_disp_sq = np.insert(cos_disp_sq, 0, 1)
    bin_edges = np.insert(bin_edges, 0, 0)
    cos_disp = np.insert(cos_disp, 0, 1)

    # Linear fit for the first two plots.
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0 : fitf], 1-cos_disp_sq[fit0 : fitf])

    # Gaussian fit for the third plot.
    b2_l = linear_fit(bin_edges[:-1]**2, *popt_linear) - (1-cos_disp)

    # Gaussian Autocorrelation Function
    popt_gauss, __ = curve_fit(gauss_function, bin_edges[:-1], b2_l)

    print("Y-intercept: ", popt_linear[-1])
    print("Amplitude, sigma")
    print("Gaussian parameters are: ", popt_gauss)
    print("FWHM: ", popt_gauss[1] * 2.35)

    # Where we display the multi fit figures.
    fig = plt.figure(num=1, figsize =(12, 12))
    plt.subplot(3, 1, 1)
    plt.title(ttl)
    plt.plot(bin_edges_sq[:-1], 1-cos_disp_sq, linestyle ="none", marker="X", label="Data Points")
    plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), label='Linear Fit', linestyle="--")
    plt.ylabel("1 - Cosdisp")
    plt.xlabel("L $^2$ (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(bin_edges[:-1], 1-cos_disp, linestyle ="none", marker="X", label="Data Points")
    plt.plot(bin_edges, linear_fit(bin_edges**2, *popt_linear), label='Linear Fit', linestyle="--")
    plt.ylabel("1 - Cosdisp")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(bin_edges[:-1], gauss_function(bin_edges[:-1],*popt_gauss), label="Gaussian Fit", linestyle="--")
    plt.plot(bin_edges[:-1], b2_l, linestyle ="none", marker="X", label="Fit and Data Difference")
    plt.ylabel("b$^2$(l)")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.show()


def data_cut(x_cen, y_cen, rad, data, show=False):
    """
    Cuts a circular region of data based on the map provided.
    """
    region = CirclePixelRegion(center=PixCoord(x=x_cen, y=y_cen), radius=rad)
    center = PixCoord(x_cen, y_cen)
    reg = CirclePixelRegion(center, rad)
    mask = reg.to_mask()
    mask = reg.to_mask(mode='center')
    dt = mask.cutout(data)
    if show:
        fig, ax = plt.subplots(figsize=(6, 6))
        region = RectanglePixelRegion(center=PixCoord(x=x_cen, y=y_cen), width=rad, height=rad)
        center = PixCoord(x_cen, y_cen)
        plt.imshow(data)
        plt.colorbar()
        region.plot(ax=ax, color='red')
        plt.show()
    return dt