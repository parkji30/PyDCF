import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits
from regions import PixCoord, CirclePixelRegion

os.chdir("/Users/a16472/Desktop/dcf_python/vcImages/")

def Stokes_Polarization(Q, U):
    """
    Calcualtes the polarization of the angle using stokes Q and U maps.
    """

    return np.arctan2(Q, U)/2


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
        x1 = (-1.0) * np.sin(angle1.reshape(1,n))
        y1 = np.cos(angle1.reshape(1,n))
        x2 = (-1.0) * np.sin(angle2.reshape(1,n))
        y2 = np.cos(angle2.reshape(1,n))

        v1 = np.array([x1, y1, np.zeros((1,n))])
        v2 = np.array([x2, y2, np.zeros((1,n))])
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
        d_ang0 = np.arctan2(np.sqrt(CdC), np.abs(vdgr)) # taking the abs here.
        return d_ang0

def cos_disp_calculations(data, ttl, ds_scale, bins, fit0=7, fitf=17):

    x, y, pix_ang, dphi = [], [], [], []
    c = 0

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not np.isnan(data[i][j]):
                x.append(i)
                y.append(j)
                pix_ang.append(data[i][j])
                # dphi.append(image_data_variance[i][j])

    x = np.array(x)
    y = np.array(y)
    ang = np.array(pix_ang)
    # dphi = np.sqrt(np.array(dphi))

    nump = len(ang) 
    w = 2.5 / 2.35 # arc seconds
    delta_r = []
    delta_phi = []
    err_dphi = []
    err_bars = []
    phi = ang

    for i in range(nump):
        delta_x_arr = x[i] - x[(i+1):(nump)]
        delta_y_arr = y[i] - y[(i+1):(nump)]
        delta_r_arr = np.sqrt(delta_x_arr**2 + delta_y_arr**2)

        sz_phi = len(delta_x_arr)
        phi_ref = np.repeat(phi[i], sz_phi)

        if len(phi_ref) > 0:
            delta_phi_arr = calc_rel_angle_crossn(phi_ref, phi[(i+1):(nump)])

        # Don't Have Errors atm.
        # err_dphi_arr = np.sqrt(dphi[i]**2 + dphi[(i+1):(nump)]**2 - 2.0 * dphi[i] * dphi[(i+1):(nump)] * np.exp((-0.25*delta_r_arr**2)/w**2))

        delta_r.append(delta_r_arr)
        delta_phi.append(delta_phi_arr)
        # err_dphi.append(err_dphi_arr)

    delta_r = np.array(delta_r)
    delta_phi = np.array(delta_phi[:-1]) # Last value is added twice. Makes sure to remove it.
    # err_dphi = np.array(err_dphi)

    delta_r = np.concatenate(delta_r).ravel() * 10 / 512 * ds_scale # CONVERT THIS TO UNITS OF PARSEC
    delta_r_squared = delta_r**2
    delta_phi = np.concatenate(delta_phi).ravel()
    # err_dphi = np.concatenate(err_dphi).ravel()
    
    nbins = bins 
    bin_edges = (np.linspace(0, nbins, 21) + 0.5) * 10 / 512 * ds_scale
    # bin_edges = np.insert(bin_edges, 0, 0)

    cos_disp, bin_edges_cos, bin_number_cos = stats.binned_statistic((delta_r), np.cos(delta_phi[:]), 'mean', bins = bin_edges)
    cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r)**2, np.cos(delta_phi[:]), 'mean', bins=bin_edges**2)
    # err_binned, bin_edges_err, bin_number_err = stats.binned_statistic(delta_r, err_dphi, 'mean', bins=bin_edges, range = (0, 1770))
    # error_brs = np.sqrt(np.sin(err_binned)**2 * err_binned**2 + (3/4) * np.cos(err_binned)**2 * err_binned**4)

    bin_edges_sq = np.insert(bin_edges_sq, 0, 0)
    cos_disp_sq = np.insert(cos_disp_sq, 0, 1)
    # err_binned = np.insert(err_binned, 0, 0)

    start = fit0
    end = fitf
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[start:end], 1-cos_disp_sq[start:end])

    cos_disp, bin_edges, bin_number_cos = stats.binned_statistic((delta_r), np.cos(delta_phi[:]), 'mean', bins = bin_edges)
    bin_edges = np.insert(bin_edges, 0, 0)
    cos_disp = np.insert(cos_disp, 0, 1)
    
    b2_l = linear_fit(bin_edges[:-1]**2, *popt_linear) - (1 - cos_disp)
    popt_gauss, __ = curve_fit(gauss_function, bin_edges[:-1], b2_l, p0=[np.max(bin_edges), np.max(bin_edges)/2]) # Gaussian Fit
    popt_gauss, __ = curve_fit(gauss_function, bin_edges_ps[:-1], b2_l, p0=popt_gauss)
    
    print("Y-intercept: ", popt_linear[-1])
    print("Amplitude, sigma")
    print("Gaussian parameters are: ", popt_gauss)
    print("FWHM: ", popt_gauss[1] * 2.35)

    return {'Bin Edges Squared': bin_edges_sq, 'Bin Edges': bin_edges, }

    fig = plt.figure(figsize =(12, 12))
    plt.subplot(3, 1, 1)
    plt.title(ttl)
    plt.plot(bin_edges_sq[:-1], 1-cos_disp_sq, linestyle ="none", marker="X", label="Data Points")
    plt.plot(bin_edges_sq[:-1], linear_fit(bin_edges_sq[:-1], *popt_linear), label='Linear Fit', linestyle="--")
    plt.ylabel("1 - Cosdisp")
    plt.xlabel("L $^2$ (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(bin_edges[:-1], 1-cos_disp, linestyle ="none", marker="X", label="Data Points")
    plt.plot(bin_edges[:-1], linear_fit(bin_edges[:-1]**2, *popt_linear), label='Linear Fit', linestyle="--")
    plt.ylabel("1 - Cosdisp")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(bin_edges[:-1], gauss_function(bin_edges[:-1], *popt_gauss), label="Gaussian Fit", linestyle="--")
    plt.plot(bin_edges[:-1], b2_l, linestyle ="none", marker="X", label="Fit and Data Difference")
    plt.ylabel("b$^2$(l)")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend()
    plt.show()