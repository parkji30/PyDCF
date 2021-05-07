import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits
from regions import PixCoord, CirclePixelRegion
from scipy.optimize import curve_fit

os.chdir("/Users/a16472/Desktop/dcf_python/testing/")

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
        d_ang0 = np.arctan2(np.sqrt(CdC), np.abs(vdgr)) #taking the abs here.
        return d_ang0


def Imshow(image, **kwargs):
    """
    Simple function to an image.
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(image)
    plt.colorbar()
    plt.show()


def pltshow(x, y):
    """
    Show relation between x and y.
    """
    plt.figure(figsize = (10, 10))
    plt.plot(x, y)
    plt.show()


def psi(Q, U):
    """
    Calculate azimuth angle.
    """

    return np.arctan2(U, Q)/2


def gauss_function(x, a, sigma):
    """
    Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2 * sigma**2))

def linear_fit(x, m, b):
    """
    Linear Fit Function.
    """
    return m*x + b


## Calculate Stokes Parameter I
# The ang returned is not scaled positively, but it really shouldn't be a problem

# I = fits.open("ccldL10B5M10n0025v2_0p65pc_PImap.fits")[0].data
# bang = fits.open("ccldL10B5M10n0025v2_0p65pc_bangmap.fits")[0].data

# Imshow(I)
bang = fits.open("/Users/a16472/Desktop/dcf_python/testing/L1M10_0.fits")[0].data
Imshow(bang)

## Looping Script that Saves
x_centers = [150, 350]
y_centers = [150, 350]

def run_image_analysis(down_sample=5, grid_space=4, rad=30):
    """
    down_sample refers to the interval frequency of sampling
    (recommend using small values for images with better resolution)

    grid_space refers to how many squares the image will be divided
    into for analysis.
    """
    for x in x_centers:
        x_cen = x
        for y in y_centers:
            y_cen = y
            # Showing how we are masking our data.
            center = PixCoord(x_cen, y_cen)
            reg = CirclePixelRegion(center, rad)
            mask = reg.to_mask() #this is how we mask the data!

            # Cut out our data here
            mask = reg.to_mask(mode='center')
            data = mask.cutout(bang)
            weighted_data = mask.multiply(bang)

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

            nump = len(ang) #-2
            w = 2.5 / 2.35 # arc seconds
            delta_r = []
            delta_phi = []
            err_dphi = []
            err_bars = []

            phi = ang
            i = 0
            while i < nump:
                delta_x_arr = x[i] - x[(i+1):(nump)]
                delta_y_arr = y[i] - y[(i+1):(nump)]
                delta_r_arr = np.sqrt(delta_x_arr**2 + delta_y_arr**2)

                sz_phi = len(delta_x_arr)
                phi_ref = np.repeat(phi[i], sz_phi)

                if len(phi_ref) > 0:
                    delta_phi_arr = calc_rel_angle_crossn(phi_ref, phi[(i+1):(nump)])

                delta_r.append(delta_r_arr)
                delta_phi.append(delta_phi_arr)

                i += down_sample # DOWNSAMPLING BY GIVEN VALUE.

            # Change to np array and Flatten
            delta_r = np.array(delta_r)
            delta_phi = np.array(delta_phi[:-1]) # last value is added twice for some reason.

            delta_r = np.concatenate(delta_r).ravel() * 10 / 512 # CONVERT THIS TO UNITS OF PARSEC
            delta_r_squared = delta_r**2
            delta_phi = np.concatenate(delta_phi).ravel()

            # FITTING CONDITIONS
            nbins = 21
            start = 10
            end = 18
            bin_edges = (np.linspace(0, nbins, nbins+1) + 0.5) * 10 / 512

            # BIN STATISTICS
            cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r)**2, \
            np.cos(delta_phi[:]),'mean', bins=bin_edges**2, range=(0, 100))
            cos_disp, bin_edges_ps, bin_number_cos = stats.binned_statistic((delta_r), np.cos(delta_phi[:]), \
            'mean', bins = bin_edges, range = (0, 10))

            bin_edges_sq = np.insert(bin_edges_sq, 0, 0)
            cos_disp_sq = np.insert(cos_disp_sq, 0, 1)
            bin_edges_ps = np.insert(bin_edges_ps, 0, 0)
            cos_disp = np.insert(cos_disp, 0, 1)

            # CURVE FITTING
            popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[start:end], 1-cos_disp_sq[start:end])
            b2_l = linear_fit(bin_edges_ps[:-1]**2, *popt_linear) - (1-cos_disp)
            popt_gauss, __ = curve_fit(gauss_function, bin_edges_ps[:-1], b2_l) #Gaussian Fit
            #popt_gauss, __ = curve_fit(gauss_function, bin_edges_am[:-1], b2_l, p0=popt_gauss)

            # FIGURES
            fig = plt.figure(figsize = (12, 8))
            plt.subplot(3, 1, 1)
            plt.title("ccldL10B5M10n0025v2_bang.fits")
            plt.plot(bin_edges_sq[:-1], 1-cos_disp_sq, linestyle ="none", marker="X", label="Data Points")
            plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), label='Linear Fit', linestyle="--")
            plt.ylabel("1 - Cosdisp")
            plt.xlabel("L $^2$ (Parsecs)")
            plt.legend()
            plt.subplot(3, 1, 2)
            plt.plot(bin_edges_ps[:-1], 1-cos_disp, linestyle ="none", marker="X", label="Data Points")
            plt.plot(bin_edges_ps, linear_fit(bin_edges_ps**2, *popt_linear), label='Linear Fit', \
            linestyle="--")
            plt.ylabel("1 - Cosdisp")
            plt.xlabel("L (Parsecs)")
            plt.legend()
            plt.subplot(3, 1, 3)
            plt.plot(bin_edges_ps[:-1], gauss_function(bin_edges_ps[:-1],*popt_gauss), label="Gaussian Fit", \
            linestyle="--")
            plt.plot(bin_edges_ps[:-1], b2_l, linestyle ="none", marker="X", label="Fit and Data Difference")
            plt.ylabel("b$^2$(l)")
            plt.xlabel("L (Parsecs)")
            plt.legend()
            plt.show()

            plt.savefig(str(x_cen) + ' ' + str(y_cen) + "ccldL10B5M10n0025v2_bang.png")

## Image analysis
run_image_analysis(down_sample=5, rad=75)


## Selecting Region or regions
x_cen = 360
y_cen = 350
rad = 30

# Showing where we are taking the cut.
fig, ax = plt.subplots(figsize=(6,6))
region = CirclePixelRegion(center=PixCoord(x=x_cen, y=y_cen), radius=rad)
ax.imshow(bang)
region.plot(ax=ax, color='red')
plt.show()

# Showing how we are masking our data.
center = PixCoord(x_cen, y_cen)
reg = CirclePixelRegion(center, rad)
mask = reg.to_mask() #this is how we mask the data!

# Cut out our data here
mask = reg.to_mask(mode='center')
data = mask.cutout(bang)
weighted_data = mask.multiply(bang)

plt.figure()
plt.title("Regional Data to Cut")
plt.imshow(data, origin='lower')
plt.colorbar()
plt.show()