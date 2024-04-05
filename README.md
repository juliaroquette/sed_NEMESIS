This repository focuses on tools related to the production of SEDs under the scope of the NEMESIS project. 

**Last Update**: @juliaroquette  - 5 April 2024

In this current version, one can import and use it by doing:

```python
import sys
sys.path.append('PAT/TO/THE/PACKAGE/LOCATION')  
```

# `flux2mag_transformations` module focus on transformations between magnitudes and fluxes.

$m=-2.5\log_{10}{\big(\frac{f}{f'}\big)},$
where $f'$ is the flux zeropoint. 

$\delta m=\frac{2.5}{\ln{10}}\frac{\delta_f}{f}$

$f = f'\times10^\frac{-m}{2.5}$

$m = m' + Ax$

$f = f'\times10^\frac{m' + Ax}{2.5} = f'\times10^\frac{m'}{2.5}10^\frac{Ax}{2.5}$



from flux2mag_transformations import SVOFilterProfileService

