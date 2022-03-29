# PyDCF (Davis Chandrasekhar Fermi Method In Python)

[![GitHub releases](https://img.shields.io/github/release/greenbone/PROJECT.svg)](https://github.com/parkji30/PyDCF/releases)
[![PyPI release](https://img.shields.io/pypi/v/PROJECT.svg)](https://pypi.org/project/PyDCF/)

## Introduction

Star formation is one of nature's many mysteries that has no clear explanation among astrophysicists and astronomers. Two fundamental forces of nature- Electromagnetism and Gravity, are believed to play a pivotal role during this process but appear to be in opposition to each other based on our current scientific theories.

The Davis Chandrasekhar Fermi method is an proposed theory that uses the polarization of light to calculate the large scale magnetic field strength of the interstellar medium. Through this method, we can analyze polarization maps in order to have a stronger understanding of the role that Magnetic fields play in the Star Formation Process.

Since it's initial proposition, the method has gone through several modificiations.

This method is based off the paper written by Houde et al. (2009)- https://arxiv.org/pdf/0909.5227.pdf.

## Installation

```python
pip3 install PyDCF
```

## Tutorial
First, we need to load the data. Go to the main file, this block of code loads the polarization data, velocity dispersion and mean density maps.

Note- you need astropy or some version of fits opener installed.

```python
data = fits.open("L1M10_0.1.fits")[0].data
velocity = fits.open("L1M10_sigmav_0.1.fits")[0].data
density = fits.open("L1M10_meanrho_0.1.fits")[0].data
```

Second, since our polarization map is 

```python
y_cen = (280)
x_cen = (140)
rad = 60

# Taking a smaller region from the entire map.
data_pol_region = data_cut(x_cen, y_cen, rad, data, show=True)
data_v_region = data_cut(x_cen, y_cen, rad, velocity, show=False)
data_rho_region = data_cut(x_cen, y_cen, rad, density, show=False)
```

