"""miri-mrs-clean.py

Identify and replace bad pixels.  Works alright for bright comet continuum.

Uses the wavelength axis to find bad pixels.
  - Slices near the first and last wavelength index are more difficult to clean,
    and artifacts may remain.

Uses a spatial median filter to replace them.

This can produce artifacts, so the data quality cube is updated to include the
bad pixels (i.e., they remain marked as bad data with dq = 513... is that OK?).

This code is known to remove line emission.  Use at your own risk.  Suggestions
for improvement are welcome.


v20230404, MSP Kelley

"""

import os
import warnings
from glob import glob
import numpy as np
import scipy.ndimage as nd
from astropy.io import fits

warnings.warn("This code is known to remove line emission.  Use at your own risk.")

force_update = False

files = sorted(
    [fn for fn in glob("processed/jw*epoch*/*_s3d.fits") if "cleaned" not in fn]
)
for fn in files:
    outf = fn.replace("_s3d", "-cleaned_s3d")

    istat = os.stat(fn)
    if os.path.exists(outf) and not force_update:
        ostat = os.stat(outf)
        if ostat.st_mtime > istat.st_mtime:
            # file is up to date
            print(fn, "- skipped")
            continue
    print(fn, "- processing")

    hdu = fits.open(fn)

    cube = hdu["SCI"].data
    cleaned = cube.copy()
    dq = hdu["DQ"].data
    for i in range(cube.shape[0]):
        im = cube[i].copy()

        # at the start and end of the cube, fewer slices should be inspected
        if i < 4:
            window = np.s_[: (i + 1) * 2]
        elif i > cube.shape[0] - 5:
            window = np.s_[-(cube.shape[0] - i - 1) * 2 - 1 :]
        else:
            window = np.s_[i - 4 : i + 5]

        # median nearby wavelengths and compare to this wavelength
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            im_wmed = np.nanmedian(cube[window], 0)
        det = (im - im_wmed) / im_wmed

        # border varies by wavelength, so identify it from the nan pixels within
        # the image, this won't be perfect when masked features touch the image
        # edge
        features, _ = nd.label(np.isnan(im))
        border = features == features[0, 0]

        # strong outliers and NaNs are masked, but preserve the border
        mask = (det > 0.5) + np.isnan(im) * ~border
        im[mask] = np.nan

        # use spatial filter to replace bad data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            im_smed = nd.generic_filter(im, np.nanmedian, 5)
        im[mask] = im_smed[mask]

        cleaned[i] = im
        dq[i][mask] = 513  # this value is just a guess

    hdu["SCI"].data = cleaned
    hdu["SCI"].header.add_history("artifacts cleaned with clean.py")
    hdu["DQ"].data = dq
    hdu.writeto(outf, overwrite=True)
