#include <stdio.h>
#include <math.h>


int linear_fit(float x, float m, float b){
    return (x * m + b);
}

int gaussian_function(float x, float a, float sigma){
    return a * exp(pow(-x, 2) / 2*pow(sigma, 2));
}
    
int total_gauss_function(float x,  float a, float W, float delta){
    return a * exp(pow(-x, 2) / (2*pow(delta, 2) + 2*pow(W, 2)));
}

int turbulent_cells(float delta, float cloud_dep, float beam_res){
    return (pow(delta, 2) + 2*pow(beam_res, 2) * cloud_dep / sqrt(2*M_1_PI) *pow(delta, 3));
}


int single_fit(data_pack, ttl, fit0, fitf, beam_res, show) {
 
    char ttl[20];
    int fit0;
    int fitf;
    float beam_res;
    bool show = true;
  
    cos_disp = data_pack[0]
    bin_edges_norm = data_pack[1]
    cos_disp_sq = data_pack[2]
    bin_edges_sq = data_pack[3]
    W = beam_res / 2.35
    
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0:fitf], 1-cos_disp_sq[fit0:fitf]) # Linear Fit
    b2_l = linear_fit(bin_edges_norm**2, *popt_linear) - (1-cos_disp) # Linear Fit Squared
    popt_gauss, __ = curve_fit(gauss_function, bin_edges_norm, b2_l) # Gaussian Fit
    
    print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", popt_linear[-1])
    print('[ Amplitude  Sigma ]')
    print("Gaussian parameters are: ", popt_gauss)
    print("FWHM: ", popt_gauss[1] * 2.35)
    
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
    plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, a=popt_gauss[0], sigma=W), label='Gaussian Beam Contribution', color='r', linestyle='--')
    plt.plot(bin_edges_norm, total_gauss_function(bin_edges_norm, popt_gauss[0], W, analytic_turb_cof), linestyle ="--", marker=".", label='Analytic Turbulent + Beam', color='g')
    plt.ylabel("b$^2$(l)")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend(loc=1)
    if show:
        plt.show()
    return popt_linear[-1], analytic_turb_cof
}

def multi_fit(data_pack, ttl, fit0, fitf, beam_res, show=True):
    """
    """
    cos_disp = data_pack[0]
    bin_edges_norm = data_pack[1]
    cos_disp_sq = data_pack[2]
    bin_edges_sq = data_pack[3]
    
    W = beam_res / 2.35
    
    popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0:fitf], 1-cos_disp_sq[fit0:fitf]) # Linear Fit
    b2_l = linear_fit(bin_edges_norm**2, *popt_linear) - (1-cos_disp) # Linear Fit Squared
    popt_gauss, __ = curve_fit(gauss_function, bin_edges_norm, b2_l) # Gaussian Fit
    
    print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", popt_linear[-1])
    print('[ Amplitude  Sigma ]')
    print("Gaussian parameters are: ", popt_gauss)
    print("FWHM: ", popt_gauss[1] * 2.35)
    
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
    plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, a=popt_gauss[0], sigma=W), label='Gaussian Beam Contribution', color='r', linestyle='--')
    plt.plot(bin_edges_norm, total_gauss_function(bin_edges_norm, popt_gauss[0], W, analytic_turb_cof), linestyle ="--", marker=".", label='Analytic Turbulent + Beam', color='g')
    plt.ylabel("b$^2$(l)")
    plt.xlabel("L (Parsecs)", fontsize=11.5)
    plt.legend(loc=1)
    if show:
        plt.show()
         
    
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