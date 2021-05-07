import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
from astropy.io import fits
from regions import PixCoord, CirclePixelRegion

os.chdir("/Users/a16472/Desktop/dcf_python/vcImages/")

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


def plt_show(x, y):
    """
    Show relation between x and y.
    """
    plt.figure(figsize = (10, 10))
    plt.plot(x, y)
    plt.show()


def hist_show(array, title='', xlabel='Value'):
    """
    Histogram of array
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

## Calculate Stokes Parameter I
# The ang returned is not scaled positively, but it really shouldn't be a problem
masking = np.load("true_masking.npy")
image_data = fits.open("VelaC_500_intermediate_regrid_30as_pix_ang.fits")[0].data
image_data_var = fits.open("VelaC_500_intermediate_regrid_30as_pix_var_ang.fits")[0].data * masking
pol_int = fits.open("VelaC_500_intermediate_regrid_30as_pix_I.fits")[0].data * masking

Imshow(image_data)
Imshow(image_data_var)
Imshow(pol_int, vmin=np.mean(pol_int)-np.std(pol_int)*3, vmax=np.mean(pol_int)+np.std(pol_int)*3)


## Selecting Region
x_cen = 60
y_cen = 60
rad = 34

# Showing where we are taking the cut.
fig, ax = plt.subplots(figsize=(6,6))
region = CirclePixelRegion(center=PixCoord(x=x_cen, y=y_cen), radius=rad)
ax.imshow(image_data)
region.plot(ax=ax, color='red')
plt.show()

# Showing how we are masking our data.
center = PixCoord(x_cen, y_cen)
reg = CirclePixelRegion(center, rad)
mask = reg.to_mask() #this is how we mask the data!

# Cut out our data here
mask = reg.to_mask(mode='center')
data = mask.cutout(image_data)
data_var = mask.cutout(image_data_var)
weighted_data = mask.multiply(image_data)

Imshow(data)
# hist_show(data.flatten(), title='Region Histogram ' + str(x_cen) + str(y_cen), \
    # xlabel='Angular Dispersion Value')
report_stats(data.flatten())

## Getting coord, angle pairs.. (Need to implement variance pair)
x, y, pix_ang, dphi = [], [], [], []

c = 0

for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        if not np.isnan(data[i][j]):
            x.append(i)
            y.append(j)
            pix_ang.append(data[i][j])
            dphi.append(data_var[i][j])

x = np.array(x)
y = np.array(y)
ang = np.array(pix_ang)
dphi = np.sqrt(np.array(dphi))


## Calculations
nump = len(ang) #-2
# pixel_scale = 0.00833333308499
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
    err_dphi_arr = np.sqrt(dphi[i]**2 + dphi[(i+1):(nump)]**2 - 2.0*dphi[i]*dphi[(i+1):(nump)] * np.exp((-0.25*delta_r_arr**2)/w**2))

    delta_r.append(delta_r_arr)
    delta_phi.append(delta_phi_arr)
    err_dphi.append(err_dphi_arr)


# Change to np array and Flatten
delta_r = np.array(delta_r)
delta_phi = np.array(delta_phi[:-1]) # last value is added twice for some reason.
err_dphi = np.array(err_dphi)

delta_r = np.concatenate(delta_r).ravel() * 30 / 60 # CONVERT THIS TO UNITS OF ARCSEC
delta_r_squared = delta_r**2
delta_phi = np.concatenate(delta_phi).ravel()
err_dphi = np.concatenate(err_dphi).ravel()


## Binned Statistics
# My calculated variance
nbins = 21

bin_edges_as = (np.linspace(0, 57, 58) + 0.5) * 30 / 60
bin_edges = (np.linspace(0, nbins, nbins+1) + 0.5) * 30 / 60

cos_disp, bin_edges_cos, bin_number_cos = stats.binned_statistic((delta_r), np.cos(delta_phi[:]), 'mean', bins = bin_edges, range = (0, 300))
cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic(delta_r**2, np.cos(delta_phi[:]), 'mean', bins=bin_edges**2, range=(0, 100))
err_binned, bin_edges_err, bin_number_err = stats.binned_statistic(delta_r, err_dphi**2, 'mean', bins=bin_edges, range = (0, 1770))
error_brs = np.sqrt(np.sin(err_binned)**2 * err_binned**2 + (3/4) * np.cos(err_binned)**2 * err_binned**4)

# Insert 0s at the end
bin_edges_sq = np.insert(bin_edges_sq, 0, 0)
cos_disp_sq = np.insert(cos_disp_sq, 0, 1)
err_binned = np.insert(err_binned, 0, 0)


## Curve Fitting (Arcmin^2)
from scipy.optimize import curve_fit

def linear_fit(x, m, b):
    return m*x + b

start = 13
end = 18

popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[start:end], 1-(cos_disp_sq/(1-0.5*err_binned))[start:end])

plt.figure(figsize = (12, 8))
plt.title("Linear Fit Arcmin $^2$")
plt.ylabel("1 - Cosdisp")
plt.xlabel("L $^2$ (Arcmin)")
plt.plot(bin_edges_sq[:-1], 1-(cos_disp_sq/(1-0.5*err_binned)), linestyle ="none", marker="X", label="Data Points")
plt.plot(bin_edges_sq[:-1], linear_fit(bin_edges_sq[:-1], *popt_linear), label='Linear Fit', linestyle="--")
plt.legend()
plt.show()

print("Calculated Y-Intercept is: ", popt_linear[1])

## Curve Fitting (Arcmin)
nbins = 21

cos_disp, bin_edges_am, bin_number_cos = stats.binned_statistic((delta_r), np.cos(delta_phi[:]), 'mean', bins = bin_edges, range = (0, 10))

bin_edges_am = np.insert(bin_edges_am, 0, 0)
cos_disp = np.insert(cos_disp, 0, 1)

plt.figure(figsize = (12, 8))
plt.title("Linear Fit Arcmin")
plt.ylabel("1 - Cosdisp")
plt.xlabel("L (Arcmin)")
plt.plot(bin_edges_am[:-1], 1-cos_disp/(1-0.5*err_binned), linestyle ="none", marker="X", label="Data Points")
plt.plot(bin_edges_am[:-1], linear_fit(bin_edges_am[:-1]**2, *popt_linear), label='Linear Fit', linestyle="--")
plt.legend()
plt.show()


## Curve Fitting (Gaussian Fit)
def gauss_function(x, a, sigma):
    """
    Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2 * sigma**2))

def houde_fitting_function(distance, b_strength, delta, effec_depth):
    """

    b_strength: turbulent to large scale magnetic field ratio
    delta: delta
    a:

    all you need to fit is a line (slope and intercept) which we can do
    easily.
    """
    W = 2.5 / 2.35 #arminutes
    # effec_depth = 6.5 #arc minutes
    return  np.sqrt(2*np.pi) * b_strength**2 * (delta**3 / ((delta**2 + 2*W**2) * effec_depth)) \
            * (np.exp(-distance**2 / (2*(delta**2 + 2*W**2))))


# func = houde_fitting_function
# b2_l = linear_fit(bin_edges_am[:-1]**2, *popt_linear) - (1-(cos_disp_sq/(1-0.5*err_binned)))
# popt_houde, __ = curve_fit(func, bin_edges_am[:-1], b2_l, p0=[0.28, 10, 3]) # Houde Fit (eq 53)
# popt_houde, __ = curve_fit(func, bin_edges_am[:-1], b2_l, p0=popt_houde)

popt_gauss, __ = curve_fit(gauss_function, bin_edges_am[:-1], b2_l, p0=[np.max(bin_edges_am), np.max(bin_edges_am)/2]) # Gauss Fit


plt.figure(figsize = (12, 8))
plt.title("Linear Fit and Data Difference")
plt.ylabel("b$^2$(l)")
plt.xlabel("L (Arcmins)")
plt.plot(bin_edges_am[:-1], b2_l, linestyle ="none", marker="X", label="Fit and Data Difference")
# plt.plot(bin_edges_am[:-1], houde_fitting_function(bin_edges_am[:-1], *popt_houde), label="Houde Fit", linestyle="--")
plt.plot(bin_edges_am[:-1], gauss_function(bin_edges_am[:-1],*popt_gauss), label="Gaussian Fit", linestyle="--")
plt.legend()
plt.show()

print("Amplitude, sigma")
print("Gaussian fit parameters are: ", popt_gauss)
print("FWHM: ", popt_gauss[1] *2.35)

# print("\nB Ratio, Delta, Effective Depth")
# print("Houde fit parameters are: ", popt_houde)


## All Three plots in one.

# Need to add a (0, 0) term to the bins:

fig = plt.figure()

plt.subplot(3, 1, 1)
plt.title("Vela C Region Analysis")
plt.plot(bin_edges_sq[:-1], 1-cos_disp_sq, linestyle ="none", marker="X", label="Data Points")
plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), label='Linear Fit', linestyle="--")
plt.ylabel("1 - Cosdisp")
plt.xlabel("L $^2$ (Arcmin)")
plt.legend()
plt.subplot(3, 1, 2)
plt.plot(bin_edges_am[:-1], 1-cos_disp, linestyle ="none", marker="X", label="Data Points")
plt.plot(bin_edges_am, linear_fit(bin_edges_am**2, *popt_linear), label='Linear Fit', linestyle="--")
plt.ylabel("1 - Cosdisp")
plt.xlabel("L (Arcmin)")
plt.legend()
plt.subplot(3, 1, 3)
plt.plot(bin_edges_am[:-1], gauss_function(bin_edges_am[:-1],*popt_gauss), label="Gaussian Fit", linestyle="--")
plt.plot(bin_edges_am[:-1], b2_l, linestyle ="none", marker="X", label="Fit and Data Difference")
plt.ylabel("b$^2$(l)")
plt.xlabel("L (Arcmin)")
plt.legend()
plt.show()

