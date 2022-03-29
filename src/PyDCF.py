import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from scipy.optimize import curve_fit
from regions import PixCoord, RectanglePixelRegion


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
        self.popt_guass = []

    def calculate_angular_dispersions(self, edge_length=0):
        """
        Calculate the angular dispersions for the L and L^2 space
        using the specified edge_length
        """

        if edge_length == 0:
            edge_length = 5 * self.beam_resolution
            print("Edge length set to 0, using 5 times the beam resolution instead...")

        self.angular_dispersion_analysis = \
                            self.angular_dispersion_calculation(self.polarization_data,
                                                            edge_length,
                                                            self.beam_resolution,
                                                            self.pixel_scale)

    def HH09_fit(self, fit0, fitf, cloud_depth):
        """
        The famous modified dcf (HH09) fitting method based off the paper
        Houde et al. 2009

        @type fit0: Int (Starting point for linear fit)
        @type fitf: Int (Finishing point for linear fit)
        @type cloud_depth: float (Thickness of molecular cloud)
        @rtype: None
        """

        def linear_fit(x, m, b):
            """
            Linear Fit Function.
            """
            return m * x + b


        def gauss_function(x, a, sigma):
            """
            Gaussian fit function.
            """
            return a * np.exp(-(x)**2 / (2 * sigma**2))


        def total_gauss_function(x, a, W, delta):
            """
            Turbulent Gaussian fit function.

            This is a simplified version of equation 21 in houde et al. 2009
            """
            return a * np.exp(-(x)**2 / (2*(delta**2 + 2*W**2)))


        def turbulent_cells(delta, cloud_dep, beam_res):
            """
            Number of turbulent cells across the line of sight.

            This is the "N" parameter in Houde et al. 2009.

            @type delta: float (Turbulent Correlation Length)
            @type cloud_dep: float (Depth of the cloud)
            @type beam_res: float (resolution of the data)
            @rtype: float

            """
            return ((delta**2 + 2*(beam_res**2)) * cloud_dep) / (np.sqrt(2*np.pi) * delta**3)


        def turbulent_autocorrelation(distance, b_ratio, delta, W, ED):
            """
            b_strength: turbulent to large scale magnetic field ratio
            delta: delta

            Used to fit the gaussian for the 3rd plot.
            """
            return  np.sqrt(2*np.pi) * b_ratio**2 * (delta**3 / ((delta**2 + 2*W**2) * ED)) \
                    * (np.exp(-distance**2 / (2*(delta**2 + 2*W**2))))


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

        print("FWHM of the Beam Contribution: ", popt_gauss[1] * 2.35)
        print("Fitted Analytic Turbulent Correlation Length: ", popt_houde[1])

        N = turbulent_cells(analytic_turb_cof/2.35, cloud_depth, W)
        uncorrected_turbulent_ratio = popt_linear[-1]

        self.turbulent_cells = N
        self.turbulent_correlation_length = analytic_turb_cof
        self.uncorrected_turbulent_ratio = uncorrected_turbulent_ratio
        self.popt_gauss = popt_gauss

    def HH09_parameters(self):
        """
        Reports the turbulent correlation length (delta) along with the
        corrected turbulent ratio (B_turb/B_large) and the turbulent contribution
        that comes from the telescope beam as defined in houde et al. 2009

        @rtype: None
        """
        N = self.turbulent_cells
        analytic_turb_cof = self.turbulent_correlation_length
        uncorrected_turbulent_ratio = self.uncorrected_turbulent_ratio

        print()
        print("Y-intercept (Uncorrected Turbulent-Ordered Ratio): ", uncorrected_turbulent_ratio)
        print('[ Amplitude  Sigma ]')
        print("Gaussian parameters are: ", self.popt_gauss)
        print("Analytic Turbulent Corrleation Length: ", analytic_turb_cof)
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
        sigma_pol = stats.circstd(self.polarization_data, high=np.pi, low=0)

        return correction_factor * np.sqrt(4*np.pi*mean_density) * (velocity_dispersion / sigma_pol)


    def SkalidisDCF_calculation(self, correction_factor=1):
        """
        Bfield calculation for the Skalidis DCF.

        everything is calculated using CGS units
        """
        mean_density = np.mean(self.density_data * (2.3 * 1.67e-24))
        velocity_dispersion = np.mean(self.velocity_data * 1e5)
        sigma_pol = stats.circstd(self.polarization_data, high=np.pi, low=0)

        return correction_factor * np.sqrt(2*np.pi*mean_density) * velocity_dispersion / np.sqrt(sigma_pol)


    def correction_factors(self, true_bfield):
        """
        return's the correction factor for the magnetic field strength.

        The correction factor is defined as Estimate / true_bfield.
        """
        classical_CF = self.ClassicalDCF_calculation() / true_bfield
        skalidis_CF = self.SkalidisDCF_calculation() / true_bfield
        HH09_CF = self.HH09DCF_calculation() / true_bfield

        print("Classical Correction Factor" , str(classical_CF))
        print("Skalidis Correction Factor" , str(skalidis_CF))
        print("HH09 Correction Factor" , str(HH09_CF))

        return [classical_CF, skalidis_CF, HH09_CF]


    def angular_dispersion_calculation(self, data, edge_length, beam_resolution, pixel_scale):
        """
        This function is used to calculate the structure function- i.e. the angular
        difference between two points as a function of distance, and use binned statistics
        to calculate the average across a certain distance range.

        NOTE- the units of measurement for edge_length, beam_resolution, pixel_scale
        should be defined by the author and consistent across each parameter.

        As in, use CGS consistently throughout the analysis or SI units depending
        on your situation.

        @type data: Array
        @type edge_length: Float (Typically this should be 5 * beam_resoultion)
        @type beam_resolution: Float (Resolution of the polarization map)
        @type pixel_scale: Float (Distance between 2 points on the map)
        @rtype: List[Arrays]
        """

        def calc_rel_angle_crossn(angle1, angle2):
            """
            Computing angular difference in code can often lead the wrapping errors-
            i.e. 359 degrees is technically the same as 1 degrees but the computer does
            not full comprehend that.

            This function is used to deal with any wrapping issues that may occur between
            Arrays angle1 and angle2.

            @type angle1: Numpy Array
            @type angle2: Numpy Array
            @rtype: Numpy Array
            """
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

        delta_r = np.concatenate(delta_r).ravel() * pixel_scale
        delta_phi = np.concatenate(delta_phi).ravel()

        bin_edge = edge_length / pixel_scale

        if beam_resolution == 0:
            nbins = 21
        else:
            nbins = np.floor((edge_length / beam_resolution * 5)) # Default use 5 samples / beam.

        print(f'Structure function analysis used: {nbins} number of bins')

        bin_edges = (np.linspace(0, bin_edge, int(nbins))) * pixel_scale
        cos_disp, bin_edges_norm, bin_number_cos = stats.binned_statistic(delta_r, np.cos(delta_phi), 'mean', bins=bin_edges)
        cos_disp_sq, bin_edges_sq, bin_number_sq = stats.binned_statistic((delta_r)**2, np.cos(delta_phi), 'mean', bins=bin_edges**2)

        cos_disp = np.insert(cos_disp, 0, 1)
        cos_disp_sq = np.insert(cos_disp_sq, 0, 1)

        return [cos_disp, bin_edges_norm, cos_disp_sq, bin_edges_sq]


    def Imshow(self, map, **kwargs):
        """
        Simple function to display an image.
        """
        if map == 'polarization':
            image = self.polarization_data
        elif map == 'velocity':
            image = self.velocity_data
        elif map == 'density':
            image = self.density_data

        plt.figure(figsize=(8, 6))
        plt.imshow(image, origin='lower')
        plt.colorbar()
        plt.show()