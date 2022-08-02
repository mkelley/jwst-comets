# Licensed with the MIT License, see LICENSE.txt
# Author: Michael S. P. Kelley
import os
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

# file name: approximate comet x, y position, 0-based index from _cal files
# (these files include the reference pixels)
files = {
    "data/jw02127001001_03101_00001_nrcb1/jw02127001001_03101_00001_nrcb1_uncal.fits": [
        742,
        1084,
    ],
    "data/jw02127001001_03101_00002_nrcb1/jw02127001001_03101_00002_nrcb1_uncal.fits": [
        749,
        1090,
    ],
    "data/jw02127001001_03101_00001_nrcblong/jw02127001001_03101_00001_nrcblong_uncal.fits": [
        1413,
        1587,
    ],
    "data/jw02127001001_03101_00002_nrcblong/jw02127001001_03101_00002_nrcblong_uncal.fits": [
        1417,
        1589,
    ],
}

pixel_readtime = 10e-6
row_overhead = 12 * pixel_readtime  #  From Bryan Hilbert, 2022-07-28
# frame_overhead = 524 * pixel_readtime  #  From Bryan Hilbert, 2022-07-28
# there's also another 10e-6 overhead before starting the next frame

rows = []
for fn, (x, y) in files.items():
    with fits.open(fn) as hdu:
        h = hdu[0].header
        t = Table(hdu["int_times"].data)
        t = t[t["integration_number"] == h["expcount"]][0]

        # account for read time of the target's center pixel
        # read outs are 512 pixels wide
        if h["fastaxis"] == 1:
            ncols = x % 512
        elif h["fastaxis"] == -1:
            ncols = (2047 - x) % 512
        elif h["fastaxis"] == 2:
            ncols = y % 512
        elif h["fastaxis"] == -2:
            ncols = (2047 - y) % 512
        else:
            raise ValueError("Unknown fastaxis", h["fastaxis"])

        if h["slowaxis"] == 1:
            nrows = x
        elif h["slowaxis"] == -1:
            nrows = 2047 - x
        elif h["slowaxis"] == 2:
            nrows = y
        elif h["slowaxis"] == -2:
            nrows = 2047 - y
        else:
            raise ValueError("Unknown slowaxis", h["slowaxis"])

        # number of pixels between first pixel and target's center pixel
        npixels = nrows * 512 + ncols

        # timing offset between first pixel and target's center pixel
        dt = npixels * pixel_readtime + nrows * row_overhead
        rows.append(
            {
                "file": os.path.basename(fn),
                "filter": h["filter"],
                "x": x,
                "y": y,
                "nrows": nrows,
                "ncols": ncols,
                "npixels": npixels,
                "dt": dt,
                "utc_mid": t["int_mid_MJD_UTC"] + dt / 86400,
                "utc_mid_iso": Time(
                    t["int_mid_MJD_UTC"] + dt / 86400, format="mjd", scale="utc"
                ).iso,
                "tdb_mid": t["int_mid_BJD_TDB"] + dt / 86400,
                "tdb_mid_iso": Time(
                    t["int_mid_BJD_TDB"] + dt / 86400, format="mjd", scale="tdb"
                ).iso,
            }
        )

tab = Table(rows)
tab.meta[
    "comments"
] = """
{}
x, y positions are 0-based, measured in _cal files.
row_overhead = {} s
pixel_readtime = {} s
dt in seconds
""".format(
    Time.now().iso, pixel_readtime, row_overhead
).splitlines()
tab.write("timing.txt", format="ascii.fixed_width_two_line", overwrite=True)
tab.pprint(-1, -1)
