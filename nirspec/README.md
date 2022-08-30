# NIRSpec

## Backgrounds

The script [backgrounds/spec2.py](backgrounds/spec2.py) processes NIRSpec observations of comet 22P/Kopff from program 1252, removing sky background and vertical striping, then processing the results with the JWST stage 2 spectroscopic pipeline.

The telescope observed the comet at all requested instrument settings, then moved to the background and repeated the observation sequence.  This is not ideal as the instrument may not return to the exact grating settings.  Typically, one would rather observe both the target and background before changing instrument settings.  However, in the interest of time efficiency, the slew to the background was only done once.

The JWST pipeline's background step compares the grating positions of the target and background files.  If they do not precisely match (as is the case for 1252), the background is not subtracted.  This script bypasses that test and subtracts a background.

In addition to the background subtraction, vertical striping (1/f noise?) is also measured with a sigma-clipped median and removed.

### Requires

* JWST pipeline
* astropy
* numpy
* scipy

### Files

* [spec2.py](backgrounds/spec2.py) - The script.
* [spec2.asdf](backgrounds/spec2.asdf) - Stage 2 configuration file with the background step disabled.
* [jwst-pipeline-log.cfg](backgrounds/jwst-pipeline-log.cfg) Logging configuration to save all INFO+ messages to spec2.log.
