"""
@juliaroquette
Created  on 9b May 2023
"""
import numpy as np
# import pandas as pd
import csv
from astropy import units as u
from astropy import constants as const


class zeropoints:
    """
    Class for loading information on the zeropoints and reference wavelength
    for different photometric systems.
    Reference: http://svo2.cab.inta-csic.es/svo/theory/fps/
    @juliaroquette 11 October 2022:
    ..........
    This reads:
    """

    def __init__(self, my_dir='/Users/juliaroquette/Work/Data/SED_builder/',
                 filename='zero_points.csv'):
        with open(my_dir + 'zero_points.csv', 'r') as f:
            reader = csv.DictReader(f, skipinitialspace=True)
            for r in reader:
                vars(self)[r['filter']] = self.band(r)

    class band:
        def __init__(self, row):
            self.lambda_ref = float(row['lambda_ref'])
            self.zp = float(row['zp'])
            self.bib = row['bib']


def flux2mag(flux, zeropoint, err_flux):
    """
    Convert fluxes to magnitudes as:
    mag = -2.5*log10(flux/zeropoint)
    err_mag = 2.5*err_flux/flux/log(10)
    ___________________________________
    input:
        flux: in Jy
        err_flux: in Jy
        zeropoint: in Jy
    """
    return -2.5*np.log10(flux/zeropoint), 2.5*err_flux/np.log(10)/flux


def mag2flux(mag, zeropoint, err_mag):
    """
    Convert magnitudes to fluxes as:
    flux = 10**(-mag/2.5)*zeropoint
    err_mag = zeropoint*ln(10)*10**(-mag/2.5)*err_mag/2.5
    ________________________________________________
    input:
        mag: in mag
        err_mag: in mag
        zeropoint: in Jy
    """
    return zeropoint*10**(- mag/2.5), zeropoint*np.log(10) * \
        (10**(- mag/2.5))*err_mag/2.5


def fluxPerUnitFrequency(lam_eff, flux, cgs=True):
    """
    Convert Fluxes to nuF_nu
    ___input___
    lam_eff -> Angstrom
    flux -> Jy
    ___output___
    nuF_nu in either cgs (if cgs=True) or W/m2
    """
    # Convert lambda eff to Angstrom
    lam_eff = lam_eff * u.Angstrom
    # Convert flux to Jansky
    flux = flux * u.Jy
    # Calculate frequency
    freq = const.c / lam_eff
    # Convert flux to erg/s/cm/cm/Hz
    flux = flux.to(u.erg / (u.s * u.cm**2 * u.Hz))
    if bool(cgs):
        # Calculate nuF_nu in erg/s/cm/cm
        nuF_nu = (flux * freq).to(u.erg / u.s / u.cm**2)
    else:
        # Calculate nuF_nu in W/m/m
        nuF_nu = (flux * freq).to(u.W / u.m**2)
    # Return the nuF_nu values as a NumPy array
    return nuF_nu.value
