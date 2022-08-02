# Astrometry

**Author**: Michael S. P. Kelley

**Date**: 2022 Aug 02

## Background

For JWST GO 2127, our team observed comet C/1995 O1 (Hale-Bopp) with NIRCam to measure its brightness and astrometric position.  The astrometry was used to test and update the JPL/Horizons orbital solution, so that we could proceed with NIRSpec observations using the 3"x3" field-of-view of the integral field unit.  This document summarizes that experience in order to help others considering similar measurements.

Details on our observation sequence can be obtained via the [information page for program GO 2127](https://www.stsci.edu/jwst/science-execution/program-information.html?id=2127) at the Space Telescope Science Institute.  A brief summary follows.  At the time of the observation, comet Hale-Bopp had a non-sidereal motion of about 3.2"/hr (0.89 mas/s).  The imaging sequence spanned a time period of 2300 s, during which the comet moved 2.0".  Two integrations were taken with a small dither motion between them.  The comet was simultaneously observed with NIRCam's B1 and B5 (AKA, Blong) detectors with the F182M and F360M filters, respectively.


## Absolute Timing

To produce useful astrometry, we needed to estimate the observation time of the comet.  Uncalibrated data products (uncal files) have an extension named INT_TIMES, which lists the start, middle, and end integration times for the first pixel read out in an array.

The NIRCam detectors take about 10.5 seconds to read out.  Given that the observations do not use a shutter, this read time produces a position dependent effective exposure time.  A source on the first pixel read out has an observation time about 10.5 seconds earlier than a source located on the last pixel.

Specific to our observations, we observed with NIRCam in full frame mode.  The detectors are simultaneously read out in four 512x2048 pixel sections, where the short axis is the read out direction (i.e., the fast axis).  The delay between the first pixel and an arbitrary pixel can be estimated knowing its position in the read out section (see [NIRCam Detector Readout](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-readout)), the fast and slow axis read directions (see [JWST Pipeline Documentation: Reference File Information](https://jwst-pipeline.readthedocs.io/en/latest/jwst/references_general/references_general.html#standard-required-keywords)), the pixel read time (10 μs), and the row read time overhead (120 μs).  With the comet near pixel (230, 963) in its read out sector, it was observed 5.05 seconds after the first pixel.  We applied this offset to the integration mid-time to produce the effective observation time for the comet.

There are several sources of timing error.  Via the JWST Help Desk, the JWST Pipeline team provided information on the absolute timing of the spacecraft clock and instrument times.  The RMS of the systematic errors in the packet time-stamps (upon which the instrument timing keywords are based) is 111 msec (three 64 msec message transfer times).  In addition, the spacecraft clock drifts from ground UTC, and the Science and Operations Center has a requirement to maintain the spacecraft clock with an accuracy of 1 second.  This difference can be examined in the [JWST Engineering Database](https://mast.stsci.edu/portal_jwst/Mashup/Clients/jwstedb/jwstedb.html).  Search for the telemetry parameter FGDP_SCTA_OFFSET_D.  It appears to have units of microseconds.  For our observations, the difference was small compared to the aforementioned 111 msec uncertainty and we did not take it into account.

A script to compute the timing offset between the first pixel and the position of comet Hale-Bopp in our data is available: [timing.py](timing.py).

## Absolute Astrometry

Astrometry was measured with the resampled 2D pipeline data products (i2d files), which are astrometric distortion corrected images.  Five Gaia DR2 astrometric sources were observed in B1, and 16 in B5 (B1 covers ~1.2 square arcmin, B5 ~4.8 square arcmin).  One of the Gaia sources was a galaxy and not used in the analysis.  The observation tracked the comet and the astrometric sources were trailed.  However, the trailing was reasonably small and linear, such that good positions could still be measured.  The source positions were fit with a linear plate solution (rotation, scale, and offset).

Fitting the image scale, rather than just applying a translational offset to align sources, could be a critical step in this process.  While the telescope tracks the comet, background stars near the first pixel are at a different location than those near the last pixel.  The result is an image that is compressed/stretched and/or skewed with respect to the true scene on the image plane.  For Hale-Bopp, the comet moves 9.5 mas per frame or about one-third of a pixel on module B1.  Compare this to our ~2 mas comet centroiding accuracy, and our 9 mas astrometric solution uncertainty (best case).

## Acknowledgements
Marco Micheli led the astrometric measurements.  Thanks to Bryan Hilbert and Tyler Pauly for providing information on the NIRCam timing and clock uncertainties.  Additional thanks to the whole JWST Hale-Bopp team.

## Notes

From Bryan Hilbert (2022 July 28):
> [A function] I wrote long ago to calculate the total exposure time of a frame: https://github.com/spacetelescope/mirage/blob/master/mirage/utils/utils.py#L93
> 
> The function can also be used to calculated the exposure time for each pixel, with a little modification. The 251ms you mentioned is made up of a 12 pixel (120 microsec) overhead associated with each row (colpad in the equation), plus a 1-row (512+12 pixel, 5240 microsec) overhead at the end of the frame. There's also a further 10 microsecond overhead at the end of the frame before moving on to the next.

From Tyler Pauly (2022 July 28):
> Generic timing accuracy and uncertainty discussion from Systems Engineer Daryl Swade:
> 
> The largest uncertainty in the FITS header timing keywords comes from the drift in the spacecraft clock. The clock is periodically updated by the flight operations team during a contact to meet the S&OC requirement to maintain an accuracy of 1 second a 24-hour period.
There is an engineering parameter, FGDP_SCTA_OFFSET_D, which tracks the difference between the spacecraft clock and UTC. The telemetry parameter FGDP_SCTA_OFFSET_D is available for any given time range from the MAST Engineering Database - https://mast.stsci.edu/portal_jwst/Mashup/Clients/jwstedb/jwstedb.html.
>
> Over the period from May 25 through June 2, 2022, FGDP_SCTA_OFFSET_D has a range of -40 to 120 msec that results from a Flight Operation Team s/c clock reset. This range is a result of a manual operation and is not a systematic error. The systematic error in FGDP_SCTA_OFFSET_D is reported to be 64 msec.
> 
> On-board, the ISIM clock is used to place a time stamp on the detector readout packets. The ISIM clock is within 128 msec of the spacecraft clock.
Each detector readout is captured in data packets within the ISIM. Each data packet is time-stamped at the time of the last pixel readout contained in the packet. In the FITS header, the ENDTIME keyword has the value of the time-stamp of the last packet of the exposure. The EXPSTART time and Barycentric times are calculated from the EXPEND time.
> 
> To summarize:
> 
> The rms of the systematic errors is 111 msec, which results from three 64 msec message transfer times between the spacecraft clock and the packet time-stamp.
The difference between the spacecraft clock and ground UTC varies over time and is set manually by the flight operations team. This difference had a range of -40 to 120 msec between May 25 and June 2, 2022.
