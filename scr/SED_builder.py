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
import pandas as pd
# list of useful markers from matplotlib
useful_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'H', '+',
                  'x', 'D','d', 'P', 'X']


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
                 filename='zero_points.csv', cmap=plt.cm.nipy_spectral):
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
        #
        survey_names = []
        for name in self.bands:
            testing_string = name.split('_')[0]
            if testing_string not in set(survey_names):
                survey_names.append(testing_string)
        self.colors = {survey_names[i]: x for i, x in
                       enumerate(cmap(np.linspace(0.05, 0.95,
                                                  len(survey_names))))}
        self.markers = {survey_names[i]: useful_markers[x] for i, x in
                        enumerate(np.random.randint(0, len(useful_markers),
                                                    len(survey_names)))}
        self.colors['generic'] = 'darkgray'

    class band:
        def __init__(self, row):
            self.lambda_ref = float(row['lambda_ref'])
            self.zp = float(row['zp'])
            # self.err_zp = float(row['err_zp'])
            self.bib = row['bib']


def plotSED(zp, df, alpha=0.9, ms=150, save=False, cmap=plt.cm.nipy_spectral,
            save_dir='/Users/juliaroquette/Work/Plots/SED/', prefix='', talkative=False, **kargs):
    """
    zp is the class zeropoints loaded, for example:
        zp = zeropoints()
    df is a dataframe containing lambdas, nuF_nu and uncertainties
    """
    try:
        assert len(df) == 1
        if bool(talkative):
            print("Correct size for {0} df: {1}".format(
                df['Internal_ID'].iloc[0], len(df)))
        df.reset_index(inplace=True)
    except AssertionError:
        print("Please, make sure df includes only one source. Current size:",
              len(df))
    nu_Fnu_key = 'nu_Fnu'
    error_key = 'error'
    lim_key = 'lim'
    column = df.columns.to_list()
    # gets all column names for nuF_nu
    nu_Fnu_columns = np.array([item for item in column if (nu_Fnu_key in item)
                              and (error_key not in item)
                              and (lim_key not in item)]) 
    if 'df2' in kargs.keys():
       df2 = kargs['df2']    
    if 'df2' in kargs.keys():
        column2 = df2.columns.to_list()        
        nu_Fnu_columns2 = np.array([item for item in column2 if (nu_Fnu_key in item)
                              and (error_key not in item)
                              and (lim_key not in item)]) 
    all_lambda = []
    all_nuFnu = []
    is_ploted = []
    fig = None
    ax = None
    if 'fig' in kargs.keys(): 
        fig = kargs['fig']
    if 'ax' in kargs.keys(): 
        ax = kargs['ax']
    if 'ax' not in kargs.keys():
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 4), dpi=200)
        plt.subplots_adjust(left=0.1, bottom=0.14, right=0.6, top=0.92)
    ax.set_yscale("log")
    ax.set_xscale("log")
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)

    # increase tick width
    ax.tick_params(width=4)
    for nm in nu_Fnu_columns:
        # print(nm)
        if pd.notna(df[nm][0]):
            all_lambda.append(df[nm[:-6] + 'lambda'][0])
            all_nuFnu.append(df[nm][0])
            survey_name = nm.split('_')[0]
            if survey_name not in list(set(is_ploted)):
                ax.scatter(df[nm[:-6] + 'lambda'][0], df[nm][0],
                        color=zp.colors[survey_name], cmap=cmap,
                        marker=zp.markers[survey_name],
                        alpha=alpha, label=survey_name, s=ms)
                is_ploted.append(survey_name)
            else:
                ax.scatter(df[nm[:-6] + 'lambda'][0], df[nm][0],
                        color=zp.colors[survey_name], cmap=cmap,
                        marker=zp.markers[survey_name],
                        alpha=alpha, s=ms)
    if 'df2' in kargs.keys():
        lb = True
        for nm in nu_Fnu_columns2:
            # print(nm)
            if pd.notna(df2[nm][0]):
                all_lambda.append(df2[nm[:-6] + 'lambda'][0])
                all_nuFnu.append(df2[nm][0])
                if bool(lb):
                    ax.scatter(df2[nm[:-6] + 'lambda'][0], df2[nm][0],
                            color='k',
                            marker='D', facecolors='none', 
                            alpha=alpha - 0.2, label='Vizier-SED', s=ms+50, zorder=0) 
                    lb = False
                else:
                    ax.scatter(df2[nm[:-6] + 'lambda'][0], df2[nm][0],
                            color='k',
                            marker='D', facecolors='none',
                            alpha=alpha - 0.2, s=ms+50, zorder=0) 
    ax.plot(np.array(all_lambda)[np.argsort(all_lambda)],
            np.array(all_nuFnu)[np.argsort(all_lambda)], 'k:', 
            alpha=0.8, zorder=0, lw=3)
    ax.set_xlabel(r'Wavelength - $\mu$m', fontsize=18)
    ax.set_ylabel(r'$\nu F(\nu)$ - erg/s/$cm^2$', fontsize=18)
    ax.set_title(f'SED for {df["Internal_ID"].iloc[0]}', fontsize=18)
    ax.set_xticks([1, 2, 10, 100, 1000])
    ax.set_xticklabels([1, 2, 10, 100, 1000])
    ax.set_xlim(0.2, 1.5e3)
    ax.set_ylim(1e-16, 1e-7)
    ax.legend(ncol=3, bbox_to_anchor=(1.05, 1), loc='upper left')
    if bool(save):
        plt.savefig(save_dir + prefix + 'SED_' + str(df["Internal_ID"].iloc[0]) +
                    '.jpg')
        plt.close(fig)
    else:
        plt.show()
    if 'fig' in kargs.keys(): 
        if 'ax' in kargs.keys(): 
            return fig, ax



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

def VizierSED2mine(id_, 
                   vizier_sed='/Users/juliaroquette/Work/Data/vizier_David/SED/raw5/reduced/'):
    """
    Function for converting SEDs obtained with the Vizier-SED web service
    into the same format as my seds.
    """
    vizier_sed = pd.read_csv(vizier_sed + f'{id_:05d}.csv')
    new_sed = pd.DataFrame({'Internal_ID' : [id_]})
    for i in range(len(vizier_sed)):
        nm = ''
        for n in vizier_sed.sed_filter[i].replace(':', ' ').replace('=', ' ').replace('/', ' ').split(' '):
            if n != ' ':
                nm += n + '_'
        new_sed[ nm + 'lambda'] = vizier_sed.sed_wl[i]
        new_sed[ nm + 'nu_Fnu']= vizier_sed.sed_fd[i]
        new_sed[ nm + 'nu_Fnu_error'] = 10**vizier_sed.log_fd_err[i]
    return new_sed
    
    