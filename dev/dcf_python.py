import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits
from regions import PixCoord, CirclePixelRegion, RectanglePixelRegion
from scipy.optimize import curve_fit
from fitting_tools import *
plt.rcParams.update({'font.size': 16})


def cos_disp_calculations(data, ds_scale):
    """
    """
    x, y, pix_ang, dphi = [], [], [], []

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


def single_fit(delta_r, delta_phi, ttl, edge_length, beam_res, fit0, fitf, beam=False, show=True):
    """
    """
    pixel_scale = 10 / 512 # User defines this but well go with Athena's for now.
    bin_edge = edge_length / pixel_scale
    nbins = np.floor(edge_length / beam_res * 5) # Always want 5 samples / beam.
    W = beam_res / 2.35
    
    bin_edges = (np.linspace(0, bin_edge, int(nbins))) * pixel_scale 

    cos_disp, bin_edges_norm, bin_number_cos = stats.binned_statistic(delta_r, np.cos(delta_phi), 'mean', bins=bin_edges)
    cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r)**2, np.cos(delta_phi), 'mean', bins=bin_edges**2)

    cos_disp = np.insert(cos_disp, 0, 1)
    cos_disp_sq = np.insert(cos_disp_sq, 0, 1)
    
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0:fitf], 1-cos_disp_sq[fit0:fitf]) # Linear Fit
    b2_l = linear_fit(bin_edges_norm**2, *popt_linear) - (1 - cos_disp) # Linear Fit Squared
    popt_gauss, __ = curve_fit(gauss_function, bin_edges_norm, b2_l) # Gaussian Fit
    
    print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", popt_linear[-1])
    print('[ Amplitude  Sigma ]')
    print("Gaussian parameters are: ", popt_gauss)
    print("FWHM: ", popt_gauss[1] * 2.35)
    print("Number of Bins: ", nbins)
    
    analytic_turb_cof = np.sqrt(popt_gauss[1]**2 - 2*(W)**2)
    print("Analytic Turbulent Corrleation Length: ", analytic_turb_cof)
    
    fig = plt.figure(num=1, figsize =(12, 12))
    plt.subplot(3, 1, 1)
    plt.title("Dispersion Analysis")
    plt.plot(bin_edges_sq[fit0:fitf],  (1-cos_disp_sq)[fit0:fitf], marker='X', label='Fitting Range', color='r')
    plt.plot(bin_edges_sq, 1-cos_disp_sq, linestyle ="none", marker=".", label=ttl + " Dispersion")
    plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), linestyle="--")
    plt.ylabel(r'$<1 - COS\phi>$')
    plt.xlabel("L $^2$ (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(bin_edges_norm[fit0:fitf],  (1-cos_disp)[fit0:fitf], marker='X', label='Fitting Range', linestyle='--', color='r')
    plt.plot(bin_edges_norm, 1-cos_disp, linestyle ="none", marker=".", label = ttl + " Dispersion")
    plt.plot(bin_edges_norm, linear_fit(bin_edges_norm**2, *popt_linear), linestyle="--")
    plt.ylabel(r'$<1 - COS\phi>$')
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, *popt_gauss), linestyle="--", label='Gaussian Fit')
    plt.plot(bin_edges_norm, b2_l, linestyle ="none", marker=".", label=ttl + r' b$^2$(l)')
    if beam:
        beam_gauss = lambda x : popt_gauss[0] * np.exp(-(x)**2 / (2*(2 * W**2)))
        plt.plot(bin_edges_norm, beam_gauss(bin_edges_norm), label='Gaussian Beam Contribution', color='r', linestyle='--')
        plt.plot(bin_edges_norm, total_gauss_function(bin_edges_norm, popt_gauss[0], W, analytic_turb_cof), linestyle ="--", marker=".", label='Analytic Turbulent + Beam', color='g')
    plt.ylabel("b$^2$(l)")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend(loc=1)
    if show:
        plt.show()
    return popt_linear[-1], analytic_turb_cof
        

def turbulent_cells(delta, cloud_dep, beam_res):
    """
    Obtains the number of turbulent cells in a beam.
    """
    return (delta**2 + 2*beam_res**2) * cloud_dep / (np.sqrt(2*np.pi)*delta**3)    
  
    
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