"""
I am following the output format provided by the SVO Filter Profile 
Service via astroquery.svo_fps

References:
- SVO documentation: http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=voservice
- astroquery.svo_fps: https://astroquery.readthedocs.io/en/latest/index.html 



TO DO:
- include a function that creates a new photometric system in the format of the SVO Filter Profile Service
- include a function that downloads the transmission curve of a photometric system
"""
import numpy as np
import pandas as pd
import astropy.units as u
import warnings

class SVOFilterProfileService:
    """
    Class for loading information on the zeropoints and reference wavelength
    for different photometric systems.
    Reference: http://svo2.cab.inta-csic.es/svo/theory/fps/
    
    @juliaroquette 4 April 2024:
    
    This version of the Class reads the Photometric System information 
    directly from the SVO database using the `astroquery.svo_fps.SvoFps`
    package. 
    """
    def __init__(self,
                nemesis_fname  = '../data/nemesis_key.csv',
                svo_fname = '../data/svo_database.csv'):
        self.data = pd.read_csv(svo_fname)
        self.data.set_index('filterID', inplace=True)
        self.nemesis2svo = pd.read_csv(nemesis_fname)
        self.nemesis2svo.set_index('filter', inplace=True)
    
    def get_filterInfo(self, band_name, nemesis_name=False):
        """
        Returns the appropriate PhotometricSystem object for the given filter (band_name).

        Args:
            band_name (_type_): _description_
            nemesis_name (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if bool(nemesis_name):
            band_name = self.nemesis2svo.loc[band_name, 'svo']
            return PhotometricSystem(self.data.loc[band_name])
        else:
            return PhotometricSystem(self.data.loc[band_name])


class PhotometricSystem:
    def __init__(self, data):
        """_summary_

        Args:
            FilterProfileService  object  Reference for the filter info
                - ivo://svo/fps	if obtained with astroquey.svo_fps
                        filterID  object  Unique identifier following the format
                                          Category/Subcategory.Filter
                  WavelengthUnit  object  Units for wavelength, typically 'Angstrom'
                   WavelengthUCD  object  Unified Content Descriptor (UCD) for SVO - em.wl
                      PhotSystem  object  Photometric system name (For example: '2MASS')
                    DetectorType  object  1 if the transmission curve gives an energy counter
                                          0  for a photon counter
                            Band  object  Photometric Band name (for example: 'J')
                      Instrument  object  Instrument name (for example: 'IRAC')
                        Facility  object  Observational Facility (for example: 'Spitzer')
                ProfileReference  object  Reference for the information about the filter (link)
            CalibrationReference  object  Reference for information about the calibration (link)
                     Description  object  A note describing that filter
                        Comments  object  Some general comment on the calculations done
        
            The following parameters are calculated from the transmission curves as described in this reference:
            - https://www.ivoa.net/documents/Notes/SVOFPS/NOTE-SVOFPS-1.0.20121015.pdf
                   WavelengthRef float64  [AA] 
                  WavelengthMean float64  [AA] Mean wavelength of the filter
                   WavelengthEff float64  [AA] Effective wavelength of the filter calculated using the 
                                               Vega spectrum as reference
                   WavelengthMin float64  [AA] First wavelength blueward where the transmission is at least 1%
                   WavelengthMax float64  [AA] Last wavelength redward where the transmission is at least 1%
                        WidthEff float64  [AA]
                   WavelengthCen float64  [AA] central wavelength between the 
                                               two wavelengths used to compute 
                                               the FWHM
                 WavelengthPivot float64  [AA] pivot wavelength
                  WavelengthPeak float64  [AA] wavelength at the maximum value of transmission
                  WavelengthPhot float64  [AA] Photon distribution wavelength
                                               calculated based on the Vega spectrum
                            FWHM float64  [AA] the difference between the two wavelengths for which
                                               filter transmission is half maximum
                            Fsun float64  [erg/(A.s.cm2)] Flux of the Sun at this wavelength
                       PhotCalID  object        filter ID followed by the system name
                          MagSys  object        Calibration system name (Vega, AB, ST, etc.)
                   ZeroPointVega float64  [Jy] Zero Point in Vega System
               ZeroPointVegaUnit  object     Units for the zero point (typically Jy)
                     ZeroPointAB float64  [Jy] Zero Point in Vega System
                 ZeroPointABUnit  object     Units for the zero point (typically Jy)
                     ZeroPointST float64  [erg/(A.s.cm2)] Zero Point in Vega System
                 ZeroPointSTUnit  object     Units for the zero point (typically Jy)                                    
                            Mag0 float64  
                   ZeroPointType  object
                       AsinhSoft float64
                TrasmissionCurve  object   Link for downloading the transmission curve
        """
        self.FilterProfileService = data['FilterProfileService']
        self.filterID = data.index
        self.WavelengthUnit = data['WavelengthUnit']
        self.WavelengthUCD = data['WavelengthUCD']
        self.PhotSystem = data['PhotSystem']
        self.DetectorType = data['DetectorType']
        self.Band = data['Band']
        self.Instrument = data['Instrument']
        self.Facility = data['Facility']
        self.ProfileReference = data['ProfileReference']
        self.CalibrationReference = data['CalibrationReference']
        self.Description = data['Description']
        self.Comments = data['Comments']
        self.WavelengthRef = data['WavelengthRef']
        self.WavelengthMean = data['WavelengthMean']
        self.WavelengthEff = data['WavelengthEff']
        self.WavelengthMin = data['WavelengthMin']
        self.WavelengthMax = data['WavelengthMax']
        self.WidthEff = data['WidthEff']
        self.WavelengthCen = data['WavelengthCen']
        self.WavelengthPivot = data['WavelengthPivot']
        self.WavelengthPeak = data['WavelengthPeak']
        self.WavelengthPhot = data['WavelengthPhot']
        self.FWHM = data['FWHM']
        self.Fsun = data['Fsun']
        self.PhotCalID = data['PhotCalID']
        self.MagSys = data['MagSys']
        self.ZeroPointVega = data['ZeroPoint']
        self.ZeroPointVegaUnit = data['ZeroPointUnit']
        #
        self.ZeroPointAB = 3631.
        self.ZeroPointVegaUnit = 'Jy'
        #
        self.ZeroPointST = 3.631e-9     
        self.ZeroPointSTUnit = 'erg/(s.cm2.A)'
        #        
        self.Mag0 = data['Mag0']
        self.ZeroPointType = data['ZeroPointType']
        self.AsinhSoft = data['AsinhSoft']
        self.TrasmissionCurve = data['TrasmissionCurve']
        
    def flux2mag(self, flux, err_flux, calibration_system='Vega'):
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
        if calibration_system == "Vega":
            zeropoint = self.ZeroPointVega
        elif calibration_system == "AB":
            zeropoint = self.ZeroPointAB
        elif calibration_system == "ST":
            zeropoint = self.ZeroPointST
        else:
            raise ValueError(f"Calibration system {calibration_system} not recognized")
        return -2.5*np.log10(flux/zeropoint), 2.5*err_flux/np.log(10)/flux
    

    
    def mag2flux(self, mag, err_mag, calibration_system='Vega'):
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
        if calibration_system == "Vega":
            zeropoint = self.ZeroPointVega
            self.flux_units = self.ZeroPointVegaUnit
        elif calibration_system == "AB":
            zeropoint = self.ZeroPointAB
            self.flux_units = self.ZeroPointABUnit
        elif calibration_system == "ST":
            zeropoint = self.ZeroPointST
            if self.ZeroPointSTUnit == 'erg/(s.cm2.A)':
                raise NotImplementedError(f"Calibration system {calibration_system} not implemented yet")
                # zeropoint = (zeropoint* u.erg / u.s / u.cm**2 / u.AA).to(u.Jy).value
                # self.flux_units = 'Jy'
        else:
            raise ValueError(f"Calibration system {calibration_system} not recognized")            
        return zeropoint * 10**(- mag/2.5), zeropoint*np.log(10) *\
            (10**(- mag/2.5))*err_mag/2.5
    

