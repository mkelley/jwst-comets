# Cube cleaning

As of July 2023, the pipeline has had substantial improvements with outlier rejection, but some artifacts or background sources may still remain.  This code attempts to identify outliers by comparing the pixel-by-pixel median along the spectral axis to the spaxel being tested.  Outliers are replaced with values from a spatial-median filter, but they remain marked as bad pixels (value 513... is that ok?) in the data quality flag.

> :warning: It is recommended to compare results before and after cleaning.  Adjust parameters as needed, and do not use if real features are being affected.  **This code is known to remove line emission.**  Suggestions for improvement are welcome.

### Requires

* astropy
* numpy
* scipy

### Files

* [miri-mrs-clean.py](miri-mrs-clean.py) - The script.
