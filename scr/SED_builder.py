"""
@juliaroquette
Created  on 9b May 2023
"""
import numpy as np
# import pandas as pd
import csv
from astropy import units as u
from astropy import constants as const
from matplotlib import colors
import matplotlib.pyplot as plt

# list of useful markers from matplotlib
useful_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', '+',
                  'x', 'D', 'd', 'P', 'X']


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
        self.bands = []
        self.lambda_refs = []
        self.zps = []
        self.refs = []
        with open(my_dir + 'zero_points.csv', 'r') as f:
            reader = csv.DictReader(f, skipinitialspace=True)
            for r in reader:
                vars(self)[r['filter']] = self.band(r)
                self.bands.append(r['filter'])
                self.lambda_refs.append(float(r['lambda_ref']))
                self.zps.append(float(r['zp']))
                self.refs.append(r['bib'])
        order = np.argsort(np.array(self.lambda_refs))
        self.bands = np.array(self.bands)[order]
        self.lambda_refs = np.array(self.lambda_refs)[order]
        self.zps = np.array(self.zps)[order]
        self.refs = np.array(self.refs)[order]

    class band:
        def __init__(self, row):
            self.lambda_ref = float(row['lambda_ref'])
            self.zp = float(row['zp'])
            # self.err_zp = float(row['err_zp'])
            self.bib = row['bib']

    def get_survey_colours(self, cmap=plt.cm.turbo, vmin=0, vmax=1):
        survey_names = []
        for name in self.bands:
            testing_string = name.split('_')[0]
            if testing_string not in set(survey_names):
                survey_names.append(testing_string)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        self.colors = {survey_names[i]: x for i, x in
                       enumerate(cmap(norm(np.linspace(0, 1,
                                                       len(survey_names)))))}
        self.markers = {survey_names[i]: useful_markers[x] for i, x in
                        enumerate(np.random.randint(0, len(useful_markers) - 1,
                                  len(self.colors)))}


def plotSED(zp, df, alpha=0.6, ms=50, save=False,
            save_dir='/Users/juliaroquette/Work/Plots/SED/'):
    """
    zp is the class zeropoints loaded, for example:
        zp = zeropoints()
    df is a dataframe containing lambdas, nuF_nu and uncertainties
    """
    lambda_key = '_lambda'
    nu_Fnu_key = '_nu_Fnu'
    error_key = 'error'
    lim_key = 'lim'
    # lim_nu_Fnu_values = []
    column = df.columns.to_list()
    nu_Fnu_values = np.array([df[item].to_list() for item in column if
                              (nu_Fnu_key in item) and (error_key not in item)
                              and (lim_key not in item)]).reshape(-1)
    lambda_values =  np.array([df[item].to_list() for item in column if
                               (lambda_key in item)]).reshape(-1)
    nu_Fnu_error_values =  np.array([df[item].to_list() for item in column if
                                     (nu_Fnu_key in item) and (error_key in
                                                               item) and
                                     (lim_key not in item)]).reshape(-1)
    zp.get_survey_colours()
    colors = [zp.colors[name.split('_')[0]] for name in column if
              (lambda_key in name)]
    markers = [zp.markers[name.split('_')[0]] for name in column if
               (lambda_key in name)]
    label = [name.split('_')[0] for name in column if (lambda_key in name)]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
    ax.set_yscale("log")
    ax.set_xscale("log")
    is_ploted = []
    for lamb, nuFnu, co, mk, nm in zip(lambda_values, nu_Fnu_values,
                                       colors, markers, label):
        if nm in set(is_ploted):
            ax.scatter(lamb, nuFnu, color=co, cmap=plt.cm.turbo, marker=mk,
                       alpha=alpha, s=ms)
        else:
            ax.scatter(lamb, nuFnu, color=co, cmap=plt.cm.turbo, marker=mk,
                       alpha=alpha, label=nm, s=ms)
            is_ploted.append(nm)
    ax.plot(lambda_values[np.argsort(lambda_values)],
            nu_Fnu_values[np.argsort(lambda_values)],
            'k:', alpha=0.8, zorder=0)
    # print(lambda_values)
    # print(nu_Fnu_values)
    ax.set_xlabel(r'Wavelength - $\mu$m')
    ax.set_ylabel(r'$\nu F(\nu)$ - erg/s/$cm^2$/Hz')
    ax.set_title(f'SED for {df["Internal_ID"].iloc[0]}')
    ax.set_xticks([0.1, 1, 2, 5, 10, 20, 100, 500, 850])
    ax.set_xticklabels([0.1, 1, 2, 5, 10, 20, 100, 500, 850])
    ax.set_xlim(0.1, 1e4)
    ax.set_ylim(1e-16, 1e-5)
    ax.legend()
    if bool(save):
        plt.savefig(save_dir + 'SED_' + str(df["Internal_ID"].iloc[0]) +
                    '.png')
        plt.close(fig)
    else:
        plt.show()


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
    return zeropoint * 10**(- mag/2.5), zeropoint*np.log(10) *\
        (10**(- mag/2.5))*err_mag/2.5


def flux2nuFnu(lam_eff, flux, error, cgs=True):
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
    # Convert flux and uncertainty to Jansky
    flux = flux * u.Jy
    error = error * u.Jy
    # Calculate frequency
    freq = const.c / lam_eff
    # Convert flux and uncertainty to erg/s/cm/cm/Hz
    flux = flux.to(u.erg / (u.s * u.cm**2 * u.Hz))
    error = error.to(u.erg / (u.s * u.cm**2 * u.Hz))
    if bool(cgs):
        # Calculate nuF_nu in erg/s/cm/cm
        nuF_nu = (flux*freq).to(u.erg / u.s / u.cm**2)
        nuF_nu_error = (error*freq).to(u.erg / u.s / u.cm**2)
    else:
        # Calculate nuF_nu in W/m/m
        nuF_nu = (flux*freq).to(u.W / u.m**2)
        nuF_nu_error = (error*freq).to(u.W / u.m**2)
    # Return the nuF_nu values as a NumPy array
    return nuF_nu.value, nuF_nu_error.value
