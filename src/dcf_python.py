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

def cos_disp_calculations(data, parameters):
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
    ds_scale = parameters[0]
    edge_length = parameters[1]
    beam_res = parameters[2]

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
    
    pixel_scale = 0.02 # User defines this but well go with Athena's for now.
    bin_edge = edge_length / pixel_scale
    
    if beam_res == 0:
        nbins = 21
    else:    
        nbins = np.floor((edge_length / beam_res * 5)) # Always want 5 samples / beam.
    
    W = beam_res / 2.35
    print(f'Structure function analysis used: {nbins} number of bins')
    
    bin_edges = (np.linspace(0, bin_edge, int(nbins))) * pixel_scale 
    cos_disp, bin_edges_norm, bin_number_cos = stats.binned_statistic(delta_r, np.cos(delta_phi), 'mean', bins=bin_edges)
    cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r)**2, np.cos(delta_phi), 'mean', bins=bin_edges**2)

    cos_disp = np.insert(cos_disp, 0, 1)
    cos_disp_sq = np.insert(cos_disp_sq, 0, 1)
        
    return [cos_disp, bin_edges_norm, cos_disp_sq, bin_edges_sq]


def single_fit(data_pack, ttl, fit0, fitf, beam_res, show=True):
    """
    """
    cos_disp = data_pack[0]
    bin_edges_norm = data_pack[1]
    cos_disp_sq = data_pack[2]
    bin_edges_sq = data_pack[3]
    W = beam_res / 2.35
    
    cos_disp[1] = (cos_disp[0] + cos_disp[2]) / 2 # Interpreting data to remove NaaN values.
    cos_disp_sq[1] = (cos_disp_sq[0] + cos_disp_sq[2]) / 2 # Interpreting data to remove NaaN values.
    
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0:fitf], 1-cos_disp_sq[fit0:fitf]) # Linear Fit
    b2_l = linear_fit(bin_edges_norm**2, *popt_linear) - (1-cos_disp) # Linear Fit Squared
    popt_gauss, __ = curve_fit(gauss_function, bin_edges_norm, b2_l) # Gaussian Fit
    popt_houde, __ = curve_fit(lambda x, b_ratio, delta: turbulent_autocorrelation(x, b_ratio, delta, W, 1.51), bin_edges_norm, b2_l)
    
    analytic_turb_cof = np.sqrt((popt_gauss[1]**2 - 2*(W)**2))
    
    if show:
        fig = plt.figure(num=1, figsize =(12, 12))
        plt.subplot(3, 1, 1)
        plt.title("Dispersion Analysis")
        plt.plot(bin_edges_sq[fit0:fitf], (1-cos_disp_sq)[fit0:fitf], marker='X', label='Fitting Range', color='r')
        plt.plot(bin_edges_sq, 1-cos_disp_sq, linestyle ="none", marker=".", label=ttl + " Dispersion")
        plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), linestyle="--")
        plt.ylabel(r'1 - $<COS\phi>$')
        plt.xlabel("L $^2$ (Parsecs)", fontsize=11.5)
        plt.legend()
        plt.subplot(3, 1, 2)
        plt.plot(bin_edges_norm[fit0:fitf],  (1-cos_disp)[fit0:fitf], marker='X', label='Fitting Range', linestyle='--', color='r')
        plt.plot(bin_edges_norm, 1-cos_disp, linestyle ="none", marker=".", label = ttl + " Dispersion")
        plt.plot(bin_edges_norm, linear_fit(bin_edges_norm**2, *popt_linear), linestyle="--")
        plt.ylabel(r'1 - $<COS\phi>$')
        plt.xlabel("L (Parsecs)", fontsize=11.5)
        plt.legend()
        plt.subplot(3, 1, 3)
        plt.plot(bin_edges_norm, turbulent_autocorrelation(bin_edges_norm, *popt_houde, W, 1.51), linestyle="--", label='Houde Turb. Function')
        plt.plot(bin_edges_norm, b2_l, linestyle ="none", marker=".", label=ttl + r' b$^2$(l)')
        plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, a=popt_gauss[0], sigma=W), label='Gaussian Beam Contribution', color='r', linestyle='--')
        plt.xlabel("L (Parsecs)", fontsize=11.5)
        plt.ylabel("b$^2$(l)")
        plt.legend(loc=1)
        plt.show()
        
        N = turbulent_cells(analytic_turb_cof/2.35, 2, W)
        uncorrected_turbulent_ratio = popt_linear[-1]
        
        print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", popt_linear[-1])
        print("Y-intercept (Uncorrected Turbulent-Ordered Ratio) Houde: ", turbulent_autocorrelation(0, *popt_houde, W, 1.51))
        print('[ Amplitude  Sigma ]')
        print("Gaussian parameters are: ", popt_gauss)
        print("FWHM: ", popt_gauss[1] * 2.35)
        print("Analytic Turbulent Corrleation Length: ", analytic_turb_cof)
        print("Fitted Analytic Turbulent Correlation Length: ", popt_houde[1])
        print()
        print("Turbulent Cells:", round(N))
        print("Corrected Turbulent to Ordered Ratio")
        print("------------------------------------")
        if N > 1:
            print("Uncorrected Turbulent Ratio:", uncorrected_turbulent_ratio)
            print("Corrected Turbulent Ratio:", np.sqrt(round(N) * uncorrected_turbulent_ratio))
            print("Corrected Turbulent Ratio ^2:", round(N) * uncorrected_turbulent_ratio)
        else:
            print(np.sqrt(uncorrected_turbulent_ratio))

    N = turbulent_cells(analytic_turb_cof/2.35, 2, W)
    uncorrected_turbulent_ratio = popt_linear[-1]
    
    return [analytic_turb_cof, N, uncorrected_turbulent_ratio]

