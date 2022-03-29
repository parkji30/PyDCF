# PyMDCF (Modified Davis Chandrasekhar Fermi Method In Python)

[![GitHub releases](https://img.shields.io/github/release/greenbone/PROJECT.svg)](https://github.com/parkji30/PyDCF/releases/tag/pilot)


Star formation is one of nature's many mysteries that has no clear explanation among astrophysicists and astronomers. Two fundamental forces of nature- Electromagnetism and Gravity, are believed to play a pivotal role during this process but appear to be in opposition to each other based on our current scientific theories.

The Davis Chandrasekhar Fermi method is an proposed theory that uses the polarization of light to calculate the large scale magnetic field strength of the interstellar medium. Through this method, we can analyze polarization maps in order to have a stronger understanding of the role that Magnetic fields play in the Star Formation Process.

Since it's initial proposition, the method has gone through several modificiations, one of which, is the famous HH09 (or MDCF) variation proposed by Martin Houde et al. (2009)- https://arxiv.org/pdf/0909.5227.pdf.


# Installation

```python
pip3 install PyDCF
```

# Tutorial
First, we need to load the data. Go to the main file, this block of code loads the polarization data, velocity dispersion and mean density maps.

Note- you need astropy or some version of fits opener installed.

```python
data = fits.open("L1M10_0.1.fits")[0].data
velocity = fits.open("L1M10_sigmav_0.1.fits")[0].data
density = fits.open("L1M10_meanrho_0.1.fits")[0].data
```

Second, since our polarization map is too large for the HH09 method, so we need to snip it down and make it smaller. Let's look at a highly filamentary region.

```python
y_cen = (280)
x_cen = (140)
rad = 60

# Taking a smaller region from the entire map.
data_pol_region = data_cut(x_cen, y_cen, rad, data, show=True)
data_v_region = data_cut(x_cen, y_cen, rad, velocity, show=False)
data_rho_region = data_cut(x_cen, y_cen, rad, density, show=False)
```

Third, load up the PyDCF package and initialize it with the data.

```python
pold1 = PyDCF(polarization = data_pol_region,
              velocity = data_v_region,
              density = data_rho_region,
              beam_resolution = 0.1,
              pixel_scale = 10/512)


pold1.calculate_angular_dispersions()
```

Finally, call the fit function as shown
```python
pold1.HH09_fit(fit0 = 18, fitf = 25, cloud_depth = 1.51)
pold1.HH09_parameters()
```

You should get a pretty plot that looks something like this!

![img1](https://user-images.githubusercontent.com/28542017/160524270-76b4520f-93c2-4f4e-8b82-07a919a35346.png)

Lastly, if we want to compare the Classical, MDCF and even the Skalidis method, we can call the follow methods:
The correction factor method takes the true magnetic field strength of the simulation and returns the estimated value divided by the true value.

```python
print(str(pold1.ClassicalDCF_calculation()*1e6) + " Microgauss")
print(str(pold1.SkalidisDCF_calculation()*1e6) + " Microgauss")
print(str(pold1.HH09DCF_calculation()*1e6) + " Microgauss")
pold1.correction_factors(10)
```
