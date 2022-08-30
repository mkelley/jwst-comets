"""
Author: Michael S. P. Kelley

Background subtract and de-stripe NIRSpec IFU data, then run the stage 2
spectroscopic pipeline.

NIRSpec observations of 22P/Kopff from program 1252 observed the comet at all
requested instrument settings, then moved to the background and repeated the
observation sequence.  This is not ideal as the instrument may not return to the
exact grating settings.  Typically, one would rather observe both the target and
background before changing instrument settings.  However, in the interest of
time efficiency, the slew to the background was only done once.

When the background step is enabled, the pipeline compares the grating positions
of the target and background files. If they do not precisely match, the
background is not subtracted.  This script bypasses that test and subtracts a
background.

In addition to the background subtraction, vertical striping (1/f noise?) is
also measured with a sigma-clipped median and removed.

"""
import os
from glob import glob
import numpy as np
import scipy.ndimage as nd
from astropy.io import fits
from astropy.stats import sigma_clip
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.background.background_sub import background_sub
from jwst import datamodels
import stdatamodels
import crds

output_dir = 'processed'

# files to process, must be _rate files
input_files = glob('data/jw01252001001_03101_0000?_nrs1/*_rate.fits')

# define the area outside of the spectra for de-striping using the sflat
h = fits.getheader(input_files[0])
ref = crds.getreferences(h, reftypes=['sflat'])
spec_mask = fits.getdata(ref['sflat']) != 0

# grow the mask by a couple pixels
spec_mask = nd.binary_dilation(spec_mask, iterations=2)

# files for manual background subtraction to bypass grating position test
background_files = glob('data/jw01252002001_03101_*nrs1/*_rate.fits')

for fn in input_files:
    # output file name
    outf = fn.replace("data/", f"{output_dir}/")
    if os.path.exists(outf):
        print('skipping', outf, "(file already exists)")
        continue

    # create directories as needed
    os.system(f'mkdir -p {os.path.dirname(outf)}')

    # copied-edited code from jwst.background.background_step
    with datamodels.open(fn) as input_model:
        bkg_model, result = background_sub(
            input_model, background_files, 3.0, None)
        result.meta.cal_step.back_sub = 'COMPLETE'

        # remove vertical stripes
        im = np.ma.MaskedArray(result.data, mask=spec_mask)
        clipped = sigma_clip(im, axis=0, sigma=2.5)
        stripes = np.outer(np.ones(im.shape[0]),
                           np.ma.mean(clipped, axis=0))
        result.data -= stripes

        # This seems like the right way to add a history entry, but I don't see
        # it in the resulting FITS file.
        comment = stdatamodels.util.create_history_entry("De-striped")
        result.history.append(comment)

        # save the file
        result.save(outf)

        # save the background
        bkg_model.save(outf.replace('_rate', '_rate_combinedbackground'))

# process with the stage 2 pipeline as usual, but do not run the background step
input_files = glob(f'{output_dir}/jw01252001001_03101_0000?_nrs1/*_rate.fits')
for fn in input_files:
    out_file = fn.replace('rate', 's3d')
    if os.path.exists(out_file):
        # compare modification times, update as needed
        rate_stat = os.stat(fn)
        s3d_stat = os.stat(out_file)
        if rate_stat.st_mtime < s3d_stat.st_mtime:
            print('skipping', out_file, '(s3d is newer than rate).')
            continue

    Spec2Pipeline.call(fn, config_file='spec2.asdf',
                       output_dir=os.path.dirname(fn),
                       logcfg='jwst-pipeline-log.cfg')
