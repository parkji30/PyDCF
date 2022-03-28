import matplotlib.pyplot as plt
from dispersion_analysis import *
from scipy.stats import circstd

class PyDCF:

    def __init__(self, polarization, velocity, density, beam_resolution, pixel_scale):
        """
        Initialize a new MDCFPy module object
        """

        self.polarization_data = polarization
        self.velocity_data = velocity
        self.density_data = density

        self.beam_resolution = beam_resolution
        self.pixel_scale = pixel_scale

        self.angular_dispersion_analysis = []

        self.turbulent_correlation_length = 0
        self.turbulent_cells = 0
        self.uncorrected_turbulent_ratio = 0

    def calculate_angular_dispersions(self, edge_length=0):
        """
        Calculate the angular dispersions for the L and L^2 space
        using the specified edge_length
        """

        if edge_length == 0:
            edge_length = 5 * self.beam_resolution
            print("Edge length set to 0, using 5 times the beam resolution instead...")

        self.angular_dispersion_analysis = \
                            angular_dispersion_calculation(self.polarization_data,
                                                            edge_length,
                                                            self.beam_resolution,
                                                            self.pixel_scale)

    def HH09_fit(self, fit0, fitf, cloud_depth):
        """
        """
        cos_disp = self.angular_dispersion_analysis[0]
        bin_edges_norm = self.angular_dispersion_analysis[1]
        cos_disp_sq = self.angular_dispersion_analysis[2]
        bin_edges_sq = self.angular_dispersion_analysis[3]

        W = self.beam_resolution / 2.35 # As defined in (Houde et al. 2009)
        cos_disp[1] = (cos_disp[0] + cos_disp[2]) / 2 # Interpreting data to remove NaaN values.
        cos_disp_sq[1] = (cos_disp_sq[0] + cos_disp_sq[2]) / 2 # Interpreting data to remove NaaN values.

        popt_linear, _ = curve_fit(linear_fit,  bin_edges_sq[fit0:fitf], 1-cos_disp_sq[fit0:fitf]) # Linear Fit
        b2_l = linear_fit(bin_edges_norm**2, *popt_linear) - (1-cos_disp) # Linear Fit Squared

        popt_gauss, __ = curve_fit(gauss_function, bin_edges_norm, b2_l) # Gaussian Fit
        popt_houde, __ = curve_fit(lambda x, b_ratio, delta: turbulent_autocorrelation(x, b_ratio, delta, W, cloud_depth), bin_edges_norm, b2_l)

        analytic_turb_cof = np.sqrt((popt_gauss[1]**2 - 2*(W)**2))

        fig = plt.figure(num=1, figsize =(12, 18))
        plt.subplot(3, 1, 1)
        plt.plot(bin_edges_sq[fit0:fitf], (1-cos_disp_sq)[fit0:fitf], marker='X', label='Fit', color='r')
        plt.plot(bin_edges_sq, 1-cos_disp_sq, linestyle ="none", marker=".", label= "L")
        plt.plot(bin_edges_sq, linear_fit(bin_edges_sq, *popt_linear), linestyle="--")
        plt.ylabel(r'1 - $<COS\phi>$', fontsize=18)
        plt.xlabel("L $^2$ [Parsec]", fontsize=18)
        plt.legend()
        plt.subplot(3, 1, 2)
        plt.plot(bin_edges_norm[fit0:fitf],  (1-cos_disp)[fit0:fitf], marker='X', label='Fit', linestyle='--', color='r')
        plt.plot(bin_edges_norm, 1-cos_disp, linestyle ="none", marker=".", label = r'L^2')
        plt.plot(bin_edges_norm, linear_fit(bin_edges_norm**2, *popt_linear), linestyle="--")
        plt.ylabel(r'1 - $<COS\phi>$', fontsize=18)
        plt.xlabel("L [Parsec]", fontsize=20)
        plt.legend()
        plt.subplot(3, 1, 3)
        plt.plot(bin_edges_norm, turbulent_autocorrelation(bin_edges_norm, *popt_houde, W, cloud_depth), linestyle="--", label='HH09 Turb. Function', color='black')
        plt.plot(bin_edges_norm, b2_l, linestyle ="none", marker=".", label=r' b$^2$(l)', color='blue')
        plt.plot(bin_edges_norm, gauss_function(bin_edges_norm, a=popt_gauss[0], sigma=W), label='Gaussian Beam Contribution', color='r', linestyle='--')
        plt.xlabel("L [Parsec]", fontsize=18)
        plt.ylabel("b$^2$(l)", fontsize=18)
        plt.legend(loc=1)
        plt.show()

        N = turbulent_cells(analytic_turb_cof/2.35, cloud_depth, W)
        uncorrected_turbulent_ratio = popt_linear[-1]

        print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", popt_linear[-1])
        print("Y-intercept (Uncorrected Turbulent-Ordered Ratio) Houde: ", turbulent_autocorrelation(0, *popt_houde, W, cloud_depth))
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

        self.turbulent_cells = N
        self.turbulent_correlation_length = analytic_turb_cof
        self.uncorrected_turbulent_ratio = uncorrected_turbulent_ratio


    def HH09DCF_calculation(self):
        """
        Bfield calculation using HH09 method.

        everything is calculated using CGS units.
        """
        mean_density = np.mean(self.density_data * (2.3 * 1.67e-24))
        velocity_dispersion = np.mean(self.velocity_data * 1e5)

        b_ratio = self.turbulent_cells * self.uncorrected_turbulent_ratio
        return np.sqrt(4 * np.pi * mean_density) * velocity_dispersion * b_ratio**(-0.5)


    def ClassicalDCF_calculation(self, correction_factor=1):
        """
        Bfield calculation for the classical DCF.

        everything is calculated using CGS units
        """
        mean_density = np.mean(self.density_data * (2.3 * 1.67e-24))
        velocity_dispersion = np.mean(self.velocity_data * 1e5)
        sigma_pol = circstd(self.polarization_data, high=np.pi, low=0)

        print(sigma_pol, mean_density, velocity_dispersion)

        return correction_factor * np.sqrt(4*np.pi*mean_density) * (velocity_dispersion / sigma_pol)


    def SkalidisDCF_calculation(self, correction_factor=1):
        """
        Bfield calculation for the Skalidis DCF.


        everything is calculated using CGS units
        """
        mean_density = np.mean(self.density_data * (2.3 * 1.67e-24))
        velocity_dispersion = np.mean(self.velocity_data * 1e5)
        sigma_pol = circstd(self.polarization_data, high=np.pi, low=0)

        return correction_factor * np.sqrt(2*np.pi*mean_density) * velocity_dispersion / np.sqrt(sigma_pol)


    def correction_factors(self, true_bfield):
        """
        return's the correction factor for the magnetic field strength.

        the correction factor is defined as Estimate / true_bfield.
        """

        print("Classical Correction Facotr" , self.ClassicalDCF_calculation / true_bfield)
        print("Skalidis Correction Facotr" , self.SkalidisDCF_calculation / true_bfield)
        print("HH09 Correction Facotr" , self.HH09DCF_calculation / true_bfield)

        return [self.ClassicalDCF_calculation / true_bfield, self.SkalidisDCF_calculation / true_bfield, self.HH09DCF_calculation / true_bfield ]


