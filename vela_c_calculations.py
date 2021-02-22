import numpy as np
import os
from scipy import stats
import matplotlib.pyplot as plt
os.chdir("C:/Users/16472/Desktop/Queens 2020 Fall/msc. Thesis/Vela C/dcf_python/data")


cos_disp_raw = np.loadtxt("cos_disp_raw.txt")
phi_values = np.loadtxt('cos_disp_true.txt')

x = cos_disp_raw.T[0]
y = cos_disp_raw.T[1]
ang = cos_disp_raw.T[2]
dphi = cos_disp_raw.T[3]

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
    err_dphi_arr = np.sqrt(dphi[i]**2 + dphi[(i+1):(nump)]**2 - 2.0*dphi[i]*dphi[(i+1):(nump)] * np.exp((-0.25*delta_r_arr**2)/w**2))

    delta_r.append(delta_r_arr)
    delta_phi.append(delta_phi_arr)
    err_dphi.append(err_dphi_arr)


# Change to np array and Flatten
delta_r = np.array(delta_r)
delta_phi = np.array(delta_phi)
err_dphi = np.array(err_dphi)

delta_r = np.concatenate(delta_r).ravel() *30
delta_r_squared = delta_r**2
delta_phi = np.concatenate(delta_phi).ravel()
err_dphi = np.concatenate(err_dphi).ravel()

delta_r_100k = delta_r[:100000]
delta_phi_100k = delta_phi[:100000]
err_dphi_100k = err_dphi[:100000]

# Look at true values
r = phi_values.T[0]
dp = phi_values.T[1]
er  = phi_values.T[2]


## Bin statistics (100k points)
bin_edges = (np.linspace(0, 57, 58) + 0.5) * 30

# My calculated variance
phi_variance, bin_edges, bin_number = stats.binned_statistic(delta_r_100k, delta_phi_100k, 'std', bins=bin_edges, range = (0, 1740))
cos_disp, bin_edges_err, bin_number_err = stats.binned_statistic(delta_r_100k, np.cos(delta_phi_100k), 'mean', bins=bin_edges, range = (0, 1740))
err_binned, bin_edges_err, bin_number_err = stats.binned_statistic(delta_r_100k, err_dphi_100k**2, 'mean', bins=bin_edges, range = (0, 1740))

true_phi_variance, bin_edges, bin_number = stats.binned_statistic(r, dp, 'std', bins=bin_edges, range = (0, 1740))
true_cos_disp, bin_edges_err, bin_number_err = stats.binned_statistic(r, np.cos(dp), 'mean', bins=bin_edges, range = (0, 1740))
true_err_binned, bin_edges_err, bin_number_err = stats.binned_statistic(r, er**2, 'mean', bins=bin_edges, range = (0, 1740))


## Plots! (100k points)
plt.figure(figsize = (12, 8))

# # # 1-COSDISP BINNED STATISTICS (NO ERROR CORRECTION)
# plt.hlines(1-cos_disp, bin_edges[:-1], bin_edges[1:], colors='r', lw=5, label='james cos disp')
# plt.hlines(1-true_cos_disp, bin_edges[:-1], bin_edges[1:], colors='g', lw=5, label='martin cos disp')
# plt.title("Mean Bin Statistics of PSI vs. DISTANCE")
# plt.ylabel("1 - Cosdisp")
# plt.xlabel("Distance Values (Arc Seconds)")

# # 1-COSDISP BINNED STATISTICS (WITH ERROR CORRECTION)
plt.hlines(1-(cos_disp/(1-0.5*err_binned)), bin_edges[:-1], bin_edges[1:], colors='r', lw=5, label='james cos disp (error corrected)')
plt.hlines(1-(true_cos_disp/(1-0.5*true_err_binned)), bin_edges[:-1], bin_edges[1:], colors='g', lw=5, label='martin cos disp (error corrected)')
plt.title("Mean Bin Statistics of PSI vs. DISTANCE")
plt.ylabel("1 - Cosdisp")
plt.xlabel("Distance Values (Arc Seconds)")

plt.legend()
plt.show()


## Binned Statistics
# My calculated variance
nbins = 21

bin_edges_as = (np.linspace(0, 57, 58) + 0.5) * 30
bin_edges = (np.linspace(0, nbins, nbins+1) + 0.5) * 30 / 60

# phi_variance, bin_edges, bin_number = stats.binned_statistic(delta_r, delta_phi[:-1], 'std', bins=bin_edges_as, range = (0, 1740))
cos_disp, bin_edges_cos, bin_number_cos = stats.binned_statistic((delta_r/60), np.cos(delta_phi[:-1]), 'mean', bins = bin_edges, range = (0, 300))
cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r/60)**2, np.cos(delta_phi[:-1]), 'mean', bins=bin_edges**2, range=(0, 100))
err_binned, bin_edges_err, bin_number_err = stats.binned_statistic(delta_r, err_dphi, 'mean', bins=bin_edges, range = (0, 1770))
error_brs = np.sqrt(np.sin(err_binned)**2 * err_binned**2 + (3/4) * np.cos(err_binned)**2 * err_binned**4)

plt.figure(figsize = (12, 8))
plt.ylabel("1 - Cosdisp")
# plt.xlabel("L$^2$ (Arc Minutes) ")
# plt.hlines(1-(cos_disp/(1-0.5*err_binned**2)), bin_edges_cos[:-1], bin_edges_cos[1:], colors='r', lw=5, label='james cos disp (error corrected 1 Million Points)')
# plt.hlines(1-(cos_disp_sq), bin_edges_sq[:-1], bin_edges_sq[1:], colors='r', lw=5, label='james cos disp (error corrected 1 Million Points)')
plt.hlines(1-(cos_disp), bin_edges[:-1], bin_edges[1:], colors='r', lw=5, label='james cos disp (error corrected 1 Million Points)')
plt.title("Mean Bin Statistics of PSI vs. DISTANCE")
plt.legend()
plt.show()


## Curve Fitting (Arcmin^2)
from scipy.optimize import curve_fit

def linear_fit(x, m, b):
    return m*x + b

start = 10
end = 18

popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[start:end], 1-cos_disp_sq[start:end])

plt.figure(figsize = (12, 8))
plt.title("Linear Fit Arcmin $^2$")
plt.ylabel("1 - Cosdisp")
plt.xlabel("L $^2$ (Arcminutes)")
plt.plot(bin_edges_sq[:-1], 1-cos_disp_sq, linestyle ="none", marker="X", label="Data Points")
plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), label='Linear Fit', linestyle="--")
plt.legend()
plt.show()

print("Calculated Y-Intercept is: ", popt_linear[1])


## Curve Fitting (Arcmin)
nbins = 20

bin_edges = (np.linspace(0, nbins, nbins+1) + 0.5) * 30 / 60
cos_disp, bin_edges_am, bin_number_cos = stats.binned_statistic((delta_r/60), np.cos(delta_phi[:-1]), 'mean', bins = bin_edges, range = (0, 10))

plt.figure(figsize = (12, 8))
plt.title("Linear Fit Arcmin")
plt.ylabel("1 - Cosdisp")
plt.xlabel("L (Arcminutes)")
plt.plot(bin_edges_am[:-1], 1-cos_disp, linestyle ="none", marker="X", label="Data Points")
plt.plot(bin_edges_am, linear_fit(bin_edges_am**2, *popt_linear), label='Linear Fit', linestyle="--")
plt.legend()
plt.show()


## Gaussian Fit
def gauss_function(x, a, sigma):
    """
    Gaussian fit function.
    """
    return a * np.exp(-(x)**2 / (2 * sigma**2))


def gauss_houde(x, ratio2, delta, W):
    """
    I believe this is the function Houde is using to apply the DCF method.

    [2., 1.5, W1, 50] The guess values
    """
    # p0=[10, 16, 1.06]
    return np.sqrt(2*np.pi) * ratio2 * (delta**3/(delta**2 + 2*W**2)*3.5) \
    * np.exp(-x**2 / 2*(delta**2 + 2 * W**2))


func = gauss_function

b2_l = linear_fit(bin_edges_am[:-1]**2, *popt_linear) - (1-cos_disp)
popt_gauss, __ = curve_fit(func, bin_edges_am[:-1], b2_l) #Gaussian Fit
popt_gauss, __ = curve_fit(func, bin_edges_am[:-1], b2_l, p0=popt_gauss)

plt.figure(figsize = (12, 8))
plt.title("Linear Fit and Data Difference")
plt.ylabel("b$^2$(l)")
plt.xlabel("L (Arcminutes)")
plt.plot(bin_edges_am[:-1], b2_l, linestyle ="none", marker="X", label="Fit and Data Difference")
plt.plot(bin_edges_am[:-1], func(bin_edges_am[:-1],*popt_gauss), label="Gaussian Fit", linestyle="--")
plt.plot(bin_edges_am[:-1], func(bin_edges_am[:-1], popt_gauss[0], (2.5/2.35)), label="Gaussian Fit (Beam Width)", linestyle="--")
plt.legend()
plt.show()

print("Amplitude, sigma")
print("Gaussian parameters are: ", popt_gauss)