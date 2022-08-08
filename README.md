# PyDCF (Davis Chandrasekhar Fermi Method In Python)

The Davis Chandrasekhar Fermi method is an proposed theory that uses the polarization of light to calculate the magnetic field strength of molecular clouds. Through this method, we can analyze polarization maps in order to have a stronger understanding of the role that magnetic fields play in the Star Formation Process. 

Many variations of the DCF method have been formulated over the past decades, one of which is the famous HH09 Analytical Dispersion method proposed by Martin Houde et al. (2009)- https://arxiv.org/pdf/0909.5227.pdf.

PyDCF serves as a Python implementation of the HH09 Analytical DCF method, a variation of the Classical DCF method which can correct for the overestimation effects not accounted in the Classical DCF method such as: beam smoothing, differential rotation, and bending of the magnetic field due to gravity.

A quick, in-depth guide is shown below. If anything is still confusing, please do not hestitate to reach out to me at jp7dec23@gmail.com. 

# Installation
PyDCF is available for installation through PyPI (Current version is 1.0.5).

```python
pip3 install PyDCF==1.0.5
```

# Tutorial
Import PyDCF to start things off. We're also going to need the fits module from astropy to load the data. If you're an Astronomer, I'm going to assume you already have astropy installed.

``` python
import PyDCF
from astropy.io import fits
```

Now, we load the data. Using the fits module, load the polarization data, velocity dispersion and mean density maps. 

You will have to change the file name in the code below to your respective files.

```python
data = fits.open("Polarization_Data.fits")[0].data
velocity = fits.open("Velocity_Dispersion_Data.fits")[0].data
density = fits.open("Mean_Density_Data.fits")[0].data
```

Since our polarization map is too large for the HH09 method to compute in a short time, lets choose to analyze a smaller region from it.

You're going to need to input the resolution and pixel scale of the data as well. Now we initialize PyDCF.

```python
pold1 = PyDCF(polarization = data_pol_region,
              velocity = data_v_region,
              density = data_rho_region,
              beam_resolution = 0.1,
              pixel_scale = 10/512)
```

:warning: **If you are dealing with large mapsr**: the computation can be quite expensive. The run-time complexity of the HH09 DCF method is O(n!) for reference. 

Let's cut the map into smaller regions for analysis instead. You can choose to just use the entire map but it could take a while.
```python
y_cen = (280)
x_cen = (140)
rad = 60

# Take a smaller region from the entire map.
data_pol_region = PyDCF.data_cut(x_cen, y_cen, rad, data, show=True)
data_v_region = PyDCF.data_cut(x_cen, y_cen, rad, velocity, show=False)
data_rho_region = PyDCF.data_cut(x_cen, y_cen, rad, density, show=False)
```

Update the data with the smaller regions. This replaces the previous polarization, velocity, density maps with the newer ones.
```python
PyDCF.update_data(data_pol_region, data_v_region, data_rho_region):
```

Next we calculate the angular dispersion function as defined by Eq. 6 in https://arxiv.org/pdf/0909.5227.pdf. This is necessary in order to calculate the magnetic field strength using the HH09 DCF method.

```python
pold1.calculate_angular_dispersions()
```

Finally, call the fit function as shown. This step is described in the method section of https://arxiv.org/pdf/0909.5227.pdf.
```python
pold1.HH09_fit(fit0 = 18, fitf = 25, cloud_depth = 1.51)
pold1.HH09_parameters()
```

You should get a pretty plot that looks something like this! 

![img1](https://user-images.githubusercontent.com/28542017/160524270-76b4520f-93c2-4f4e-8b82-07a919a35346.png)


Now you can calculate the magnetic field strength through using the HH09 DCF method. I've also included the Classical DCF (https://articles.adsabs.harvard.edu/pdf/1953ApJ...118..113C) and Skalidis DCF (https://arxiv.org/pdf/2010.15141.pdf) methods in this Python Package.

The DCF methods return the field strength in Gauss.

```python
print(str(pold1.ClassicalDCF_calculation()) + " Gauss")
print(str(pold1.SkalidisDCF_calculation()) + " Gauss")
print(str(pold1.HH09DCF_calculation()) + " Gauss")
```


