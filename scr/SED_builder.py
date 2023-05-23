"""
@juliaroquette
Created  on 9b May 2023
"""
import numpy as np
import pandas as pd
import csv

class zeropoints:
    """
    Class for loading information on the zeropoints and reference wavelength
    for different photometric systems.
    Reference: http://svo2.cab.inta-csic.es/svo/theory/fps/
    
    @juliaroquette 11 October 2022:
    
    This reads 
    """
    def __init__(self, my_dir='/Users/juliaroquette/Work/Data/SED_builder/', \
                filename='zero_points.csv'):
        with open(my_dir + 'zero_points.csv', 'r') as f:
            reader = csv.DictReader(f, skipinitialspace=True)        
            for r in reader: 
                vars(self)[r['filter']] = self.band(r)
    class band:
        def __init__(self, row):
            self.lambda_ref = float(row['lambda_ref'])
            self.zp = float(row['zp'])
            #self.err_zp = float(row['err_zp'])
            self.bib = row['bib']
            
def flux2mag(flux, zeropoint, err_flux):
    """
    Convert fluxes to magnitudes as:
    mag = -2.5*log10(flux/zeropoint)
    err_mag = 2.5*err_flux/flux/log(10)
    
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
    
    input:
        mag: in mag
        err_mag: in mag
        zeropoint: in Jy
    """
    return zeropoint*10**(- mag/2.5), zeropoint*np.log10(10)*(10**(- mag/2.5))*err_mag/2.5