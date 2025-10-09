# jwst-comets

Tools, tips, tricks for comet observations with JWST. Some may even be useful for asteroids.

- For Exposure Time Calculator tips and a tool for estimating comet brightness, see [estimations](estimations/).
- For Astronomer's Proposal Tool (APT) tips and a tool for visualizing moving targets in APT, see [apt-targets](apt-targets/).
- For some tips measuring astrometry with NIRCam, see [astrometry](astrometry/).
- To plot the fields of view of JWST data with respect to a moving target position, see [instrument-fov](instrument-fov/).
- For some tips on processing NIRSpec data, see [nirspec](nirspec/).
- I tried to do some cube cleaning on MIRI data with mixed results, see [cube-cleaning](cube-cleaning/)
- To find comets in the MAST Duplicate Observations tool, search with column [LEVEL1=comet](https://mast.stsci.edu/duplication/#/results?resolve=true&radius=3&radius_units=arcminutes&template_types=MIRI+Coron,MIRI+CoronCal,MIRI+ExtFlat,MIRI+Imaging,MIRI+LRS,MIRI+MRS,NIRCam+Coron,NIRCam+Engin,NIRCam+Grism+TS,NIRCam+Imaging,NIRCam+ImageIPR,NIRCam+ImagePIL,NIRCam+TS,NIRCam+WFSS,NIRISS+AMI,NIRISS+ExtCal,NIRISS+Imaging,NIRISS+IntFlat,NIRISS+SOSS,NIRISS+WFSS,NIRSpec+BOTS,NIRSpec+FSS,NIRSpec+Focus,NIRSpec+IFU,NIRSpec+Imaging,NIRSpec+IntLamp,NIRSpec+Anneal,NIRSpec+MOS,NIRSpec+MIMF,FGS+ExtCal,FGS+IntFlat,Realtime,WFSC+Commission,WFSC+Global,WFSC+Coarse,WFSC+Fine&status_types=COMPLETED,SKIPPED,WITHDRAWN,SCHEDULED,FAILED,IMPLEMENTATION,FLIGHT_READY&custom_col_sel_0=level_1&custom_col_val_0=comet&select_cols=dec_computed,dec_proper_motion,epoch,instrument,observation,optical_config,phase_coverage,phase_critical_period,program,category,ra_computed,ra_proper_motion,standard_target_name,starting_phase_window_midpoint,target_id,target_name,target_type,template,time_series_observation,visit,visit_status,level_1&search_key=079620b3b11d&useStore=false)

## References

File naming scheme:
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html>

Data products info:
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/index.html>

Science product info (e.g., FITS extensions):
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/science_products.html>

Calibration files at:
<https://jwst-crds.stsci.edu/>

Updated NIRCam calibration constants, ca. Aug 2022: [Boyer et al. 2022](https://arxiv.org/abs/2209.03348) ([data behind the figure](https://ers-stars.github.io/data/boyer2022_nircamoffsets.txt)).

Author: Michael S. P. Kelley
