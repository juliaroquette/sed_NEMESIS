This repository focuses on tools related to the production of SEDs under the scope of the NEMESIS project. 

**Last Update**: @juliaroquette  - 8 April 2024



# `flux2mag_transformations` module focus on transformations between magnitudes and fluxes.

Magnitudes are defined as $m=-2.5\log_{10}{\big(\frac{f}{f'}\big)},$ where $f'$ is the flux zero-point. Fluxes and zero-points are in Jy units and magnitudes are in mag. 

Uncertainties are propagated as:

$\delta m=\frac{2.5}{\ln{10}}\frac{\delta_f}{f}$

To revert, fluxes can be retrieved from magnitudes as $f = f'\times10^\frac{-m}{2.5}$ as  $f'\times\log_{10}(10^{- \frac{m}{2.5}})\times\frac{\delta m}{2.5}$. We are ignoring the asymmetry of the uncertainty propagation.

Note that so far $m$ are observed magnitudes and $f$ are observed fluxes. In this context, observed magnitudes can be broken down in:

$m = m' + A_x + 5\log_{10}\Big(\frac{d}{10 \mathrm{pc}}\Big)$, where #m' is the absolute magnitude, $A_x$ is the extinction contribution at the photometric band $x$ and $5\log_{10}\Big(\frac{d}{10 \mathrm{pc}}\Big)$ is the distance modulus (DM) to the source. Accordingly:

$f = f'\times10^\frac{m' + A_x + DM}{2.5} = f'\times10^\frac{m'}{2.5}10^\frac{A_x}{2.5}0^\frac{DM}{2.5}$


## Usage:

In this current version, one can import and use it by doing:

```python
import sys
sys.path.append('PAT/TO/THE/PACKAGE/LOCATION')  
```

To load the our offline version of the SVO database:

```python 
from flux2mag_transformations import SVOFilterProfileService

svo_object = SVOFilterProfileService()
```

To load the information for a specific photometric band:


```python 
zp = test.get_filterInfo('Generic/Johnson')
```

To print some info for this filter, for example, the filter ID on SVO and its effective wavelength:

```python
 print(zp.filterID , zp.WavelengthEff)
```

To use the SVO info to convert from magnitudes to fluxes:

```python
   flux, flux_error = zp.mag2flux(mag, mag_error)
```

In this example, `mag` and  `mag_error` are the magnitude and its uncertainty in the `Generic/Johnson.U` band.

If the filter name follows the NEMESIS Orion compilation nomenclature, the zero-point info can be loaded using the NEMESIS column names as long as `nemesis_name=True`. Valid NEMESIS names are listed in the `data/nemesis_key.csv` table. In this example, `Generic/Johnson.U` is called `johnson_U` within the NEMESIS catalogue. 

```python
    zp = svo.get_filterInfo('johnson_U', nemesis_name=True)
```

