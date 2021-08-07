import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from fitting_tools import *

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


def cos_disp_calculations(data, pixel_unit, W):
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
    delta_r = np.concatenate(delta_r).ravel() * pixel_unit
    delta_phi = np.concatenate(delta_phi).ravel()
    return delta_r, delta_phi


def MDCF_fit(data, pixel_scale, edge_length, beam_size, fit0, fitf, name='Region'):
    """
    """
    delta_r, delta_phi = cos_disp_calculations(data, pixel_scale, beam_size)
    bin_edge = edge_length / pixel_scale
    nbins = np.floor((edge_length / (beam_size * 2.35)) * 5 ) # Always want 5 samples / beam.
    
    print(nbins)
    W = beam_size
    bin_edges = (np.linspace(0, bin_edge, int(nbins))) * pixel_scale 

    cos_disp, bin_edges_norm, bin_number_cos = stats.binned_statistic(delta_r, np.cos(delta_phi), 'mean', bins=bin_edges)
    cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r)**2, np.cos(delta_phi), 'mean', bins=bin_edges**2)
    cos_disp = np.insert(cos_disp, 0, 1)
    cos_disp_sq = np.insert(cos_disp_sq, 0, 1)

    # Fits for subplots
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0:fitf], 1-cos_disp_sq[fit0:fitf]) # Linear Fit
    b2_l = linear_fit(bin_edges_norm**2, *popt_linear) - (1 - cos_disp) # Linear Fit Squared
    popt_gauss, __ = curve_fit(gauss_function, bin_edges_norm, b2_l) # Gaussian Fit
    
    print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", popt_linear[-1])
    print('[ Amplitude  Sigma ]')
    print("Gaussian parameters are: ", popt_gauss)
    print("FWHM: ", popt_gauss[1] * 2.35)
    print("Number of Bins: ", nbins)
    
    # Turbulent Correlation Length
    analytic_turb_cof = np.sqrt(popt_gauss[1]**2 - 2*(W)**2)
    print("Analytic Turbulent Corrleation Length: ", analytic_turb_cof)
    
    fig = plt.figure(num=1, figsize =(12, 12))
    plt.subplot(3, 1, 1)
    plt.title("Dispersion Analysis")
    plt.plot(bin_edges_sq[fit0:fitf],  (1-cos_disp_sq)[fit0:fitf], marker='X', label='Fitting Range', color='r')
    plt.plot(bin_edges_sq, 1-cos_disp_sq, linestyle ="none", marker=".", label=name + " Dispersion")
    plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), linestyle="--")
    plt.ylabel(r'$<1 - COS\phi>$')
    plt.xlabel("L $^2$ (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(bin_edges_norm[fit0:fitf],  (1-cos_disp)[fit0:fitf], marker='X', label='Fitting Range', linestyle='--', color='r')
    plt.plot(bin_edges_norm, 1-cos_disp, linestyle ="none", marker=".", label = name + " Dispersion")
    plt.plot(bin_edges_norm, linear_fit(bin_edges_norm**2, *popt_linear), linestyle="--")
    plt.ylabel(r'$<1 - COS\phi>$')
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, *popt_gauss), linestyle="--", label='Gaussian Fit')
    plt.plot(bin_edges_norm, b2_l, linestyle ="none", marker=".", label=name + r' b$^2$(l)')
    plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, a=popt_gauss[0], sigma=W), label='Gaussian Beam Contribution', color='r', linestyle='--')
    # plt.plot(bin_edges_norm, total_gauss_function(bin_edges_norm, popt_gauss[0], W, analytic_turb_cof), linestyle ="--", marker=".", label='Analytic Turbulent + Beam', color='g')
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.ylabel("b$^2$(l)")
    plt.legend(loc=1)
    return popt_linear[-1], analytic_turb_cof
        

def turbulent_cells(delta, cloud_dep, beam_res):
    """
    Obtains the number of turbulent cells in a beam.
    """
    return (delta**2 + 2*beam_res**2) * cloud_dep / (np.sqrt(2*np.pi)*delta**3)    
  
  
def staikos_DCF(field_density, sigma_v, sigma_pol):
    """
    The Modified DCF method as written by Staikos et al. 2021
    """
    return np.sqrt(2*np.pi*field_density) * sigma_v / np.sqrt(sigma_pol)


def classical_DCF(field_density, sigma_v, sigma_pol):
    """
    The classical dcf method as written by David, Chandrasekhar and Fermi.
    """
    return np.sqrt(4*np.pi*field_density) * sigma_v / sigma_pol


def modified_DCF(field_density, velocity_los, b_ratio):
    """
    The modified DCF method as written by Houde et al. (2008, 2011, 2013, 2016)
    """
    return np.sqrt(4*np.pi*field_density) * velocity_los * b_ratio**(-0.5)

