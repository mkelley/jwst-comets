# jwst-comets
Tools, tips, tricks for designing comet observations with JWST

## Tips
### Power-law profile

Cometary comae have surface brightness profiles that tend to follow a radial power-law.  The ideal coma is proportional to 1/ρ, where ρ is the distance to the nucleus.  The JWST ETC has a power-law surface brightness profile.  The central core of that profile is flat, then it follows 1/ρ thereafter.  This can be used for comets, but with caution.  In order to faithfully reproduce a cometary scene, the central core must be smaller than the pixel size of the instrument in question.

Based on the script `etc-profile-test.py` and the notebook `notes/etc-profile-notes.py`, I recommend the following:
* Generate the surface brightness of a coma averaged over an aperture of radius of 2 mas in units of mJy/arcsec^2.
* Upload this spectrum to the JWST ETC.
* Create a new source, using this spectrum as the continuum.
* Set "Do not renormalize" under the Renorm tab.
* Set the shape to "Power law":
  - normalize to per square arcsec
  - core radius to 0.001 arcsec (= 1 mas)
  - power law index to 1

There is a difference of a factor of 2 between the core radius and the aperture radius used in the coma flux estimation tool.  This difference sets the correct normalization for the ETC coma surface brightness.

## estimations
### jwst-comet-est.py

---
**Warning**

* Conversions from total visual magnitude to dust Afρ is experimental.

* Gas output is not yet implemented.

---


```
usage: jwst-comet-est.py [-h] [--aper APER] [-R R] [--unit UNIT] [-m]
                         [--dusty] [--gassy] (--afrho AFRHO | -Q Q)
                         [--ef2af EF2AF] [--Tscale TSCALE] [--h2o H2O]
                         [--co2 CO2] [--co CO] [-o O]
                         rh delta

Generate model comet spectra for JWST ETC.

positional arguments:
  rh               heliocentric distance (au)
  delta            observer-target distance (au)

optional arguments:
  -h, --help       show this help message and exit
  --aper APER      aperture radius (arcsec), or one of: ifu, fs200, fs400,
                   fs1600, lrs, mrs, nircam (default: 0.2 arcsec)
  -R R             spectral resolution (default: 10)
  --unit UNIT      output spectral flux density or surface brightness unit
                   (default: mJy/arcsec2)
  -m               indicates that --afrho or -Q is a total visual mangitude
                   that should be converted before using (default: False)
  --dusty, -d      indicates that the conversion from -m should assume a dusty
                   coma, repeated options increase dust/gas: Q/1.6|2.4|3.7,
                   Afrho*4.8|23|110 (default: 0)
  --gassy, -g      indicates that the conversion from -m should assume a gassy
                   coma, repeated options increase gas/dust: Q*1.6|2.4|3.7,
                   Afrho/4.8|23|110 (default: 0)
  --afrho AFRHO    dust mode based on this Afρ at 0 deg phase (cm) or total
                   visual magnitude (requires -m) (default: None)
  -Q Q             gas mode based on this scale factor (molecules/s) or total
                   visual magnitude (requires -m) (default: None)
  -o O             save to this file name (default: None)

dust options:
  --ef2af EF2AF    the ratio εfρ/Afρ (default: 3.5)
  --Tscale TSCALE  LTE blackbody temperature scale factor (default: 1.1)

gas options:
  --h2o H2O        relative H2O production rate (default: 1.0)
  --co2 CO2        relative CO2 production rate (default: 0.15)
  --co CO          relative CO production rate (default: 0.05)

For gas output, -m uses the Jorda et al. (2008, ACM, 8046) correlation
Q=10**(30.675 - 0.2453 * mH), where mH is heliocentric magnitude. This assumes
--h2o=1.0. --gassy and --dusty can be repeatedly used to scale the results.
```

Generate a spectrum for a comet with Afρ=100 cm at rh=1.5 au, and Δ=1.0 au, as observed by MIRI LRS:
```
$ python3 jwst-comet-est.py 1.5 1.0 --afrho=100 --aper=lrs --unit=mJy
# %ECSV 0.9
# ---
# datatype:
# - {name: wave, unit: um, datatype: float64}
# - {name: total, unit: mJy, datatype: float64}
# - {name: F_sca, unit: mJy, datatype: float64}
# - {name: F_th, unit: mJy, datatype: float64}
# meta: !!omap
# - {cmd: estimations/jwst-comet-est.py 1.5 1.0 --afrho=100 --aper=lrs --unit=mJy}
# - {rh: 1.5 AU}
# - {delta: 1.0 AU}
# - {phase: 41.40962210927086 deg}
# - {aper: "Rectangular aperture, dimensions 0.5\xD71.6 arcsec"}
# - phase function: !numpy.ndarray
#     buffer: !!binary |
#       aUZRdnJlZjgxajg9
#     dtype: float64
#     order: C
#     shape: !!python/tuple [1]
# - {ef2af: 3.5}
# - {Tscale: 1.1}
# - {Afrho: 100.0 cm}
# - {efrho: 350.0 cm}
# schema: astropy-2.0
wave total F_sca F_th
0.5 0.09409975009087415 0.09409975009087415 1.088440633433928e-41
0.55 0.110780768520771 0.110780768520771 1.93395762848993e-37
...
```
