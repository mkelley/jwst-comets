#!/usr/bin/env python3
"""Plots the S_REGION keywords for each file."""
import argparse
from copy import deepcopy
from itertools import cycle
from astropy.time import Time
from spherical_geometry.polygon import (
    SphericalPolygon,
    great_circle_arc,
    vector,
)

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord

from astroquery.jplhorizons import Horizons


parser = argparse.ArgumentParser()
parser.add_argument(
    "files",
    nargs="+",
    help="JWST data in FITS format, S_REGION keyword is required",
)
parser.add_argument(
    "--target",
    "-t",
    help="target comet or asteroid name, must be resolvable by JPL Horizons",
)
parser.add_argument(
    "--id-type",
    choices=["comet", "designation", "smallbody", "majorbody"],
    default="smallbody",
    help="target ID type (use comet to avoid matching multiple cometary orbits)",
)
parser.add_argument(
    "--no-cache",
    dest="cache",
    action="store_false",
    help="do not use cached ephemeris when possible",
)
parser.add_argument("-o", help="save plot to this file")
args = parser.parse_args()


def get_fovs(headers):
    times = {"start": [], "mid": [], "end": []}
    fovs = []
    for h in headers:
        if "S_REGION" not in h:
            raise ValueError("Expected S_REGION FITS header keyword not found.")

        fov = h["S_REGION"]
        if "POLYGON ICRS" not in fov:
            raise ValueError(
                f'Expected "POLYGON ICRS" not found in S_REGION string: {fov}'
            )

        fov = [float(x) for x in fov.split()[2:]]
        ra = fov[::2]
        dec = fov[1::2]
        fovs.append(SkyCoord(ra, dec, unit="deg"))

        times["start"].append(h["MJD-BEG"])
        times["mid"].append(h["MJD-AVG"])
        times["end"].append(h["MJD-END"])

    for k, v in times.items():
        times[k] = Time(v, format="mjd", scale="utc")

    return times, fovs


def get_ephemeris(epochs, args):
    opts = dict(cache=args.cache)

    # comets need special care: use closest orbital elements by date, do not
    # search for related fragments
    if args.id_type == "comet":
        id_type = "designation"
        opts.update(closest_apparition=True, no_fragments=True)
    else:
        id_type = args.id_type

    horizons = Horizons(
        args.target, location="@jwst", epochs=epochs.jd, id_type=id_type
    )
    eph = horizons.ephemerides(**opts)
    eph["date"] = Time(eph["datetime_jd"], format="jd")

    return eph


def quad_to_poly(ra, dec, **kwargs):
    # from catch
    p = SphericalPolygon.from_radec(ra, dec, degrees=True)

    points = p.polygons[0]._points
    _ra, _dec = [], []
    for A, B in zip(points[0:-1], points[1:]):
        length = great_circle_arc.length(A, B, degrees=True)
        if not np.isfinite(length):
            length = 2
        interpolated = great_circle_arc.interpolate(A, B, length * 4)
        lon, lat = vector.vector_to_lonlat(
            interpolated[:, 0],
            interpolated[:, 1],
            interpolated[:, 2],
            degrees=True,
        )
        for lon0, lat0, lon1, lat1 in zip(
            lon[0:-1], lat[0:-1], lon[1:], lat[1:]
        ):
            _ra.append(lon0)
            _dec.append(lat0)

    _ra.append(lon1)
    _dec.append(lat1)
    return PolyCollection([np.c_[_ra, _dec]], **kwargs)


def plot_fovs(ax, fovs, colors):
    for i, fov in enumerate(fovs):
        poly = quad_to_poly(
            fov.ra.deg, fov.dec.deg, color=next(colors), alpha=0.33
        )
        ax.add_artist(poly)
        # ax.annotate(
        #     str(i),
        #     [fov.ra.deg.mean(), fov.dec.deg.mean()],
        #     color="k",
        #     alpha=0.5,
        # )


def plot_moving_target(ax, times, args, colors):
    t = Time(np.r_[times["start"][0], times["mid"], times["end"][-1]])
    eph = get_ephemeris(t, args)
    ra = eph["RA"]
    dec = eph["DEC"]
    ax.plot(ra, dec, color="k", lw=1, label=args.target)
    ax.scatter(
        ra[1:-1],
        dec[1:-1],
        marker="x",
        c=[next(colors) for i in range(len(t) - 2)],
    )

    angles = [
        (eph["sunTargetPA"].mean() - 180) % 360,
        eph["velocityPA"].mean() - 180,
    ]
    for angle, label in zip(angles, (r"$\odot$", "v")):
        length = 0.1  # axes fraction
        a = np.radians(angle)
        dxy = length * np.array([-np.sin(a), np.cos(a)])
        xy = np.array((0.13, 0.13))
        ax.annotate(
            label,
            xy,
            xy + dxy,
            color="k",
            ha="center",
            va="center",
            fontsize=12,
            arrowprops=dict(arrowstyle="<-", shrinkB=0),
            xycoords="axes fraction",
        )

    return ra, dec


headers = sorted(
    [fits.getheader(f, ext=("sci", 1)) for f in args.files],
    key=lambda h: h["mjd-avg"],
)
times, fovs = get_fovs(headers)

fig = plt.figure(1, (8, 8), clear=True)
ax = fig.add_subplot()
colors = cycle(mpl.colors.TABLEAU_COLORS)
plot_fovs(ax, fovs, deepcopy(colors))
if args.target:
    target_ra, target_dec = plot_moving_target(
        ax, times, args, deepcopy(colors)
    )
else:
    target_ra = []
    target_dec = []

xlim = (
    np.max([fov.ra.deg.max() for fov in fovs] + list(target_ra)),
    np.min([fov.ra.deg.min() for fov in fovs] + list(target_ra)),
)
ylim = (
    np.min([fov.dec.deg.min() for fov in fovs] + list(target_dec)),
    np.max([fov.dec.deg.max() for fov in fovs] + list(target_dec)),
)
pad = 0.05 * (ylim[1] - ylim[0])
ylim = [ylim[0] - pad, ylim[1] + pad]

# adjust RA given Decl
ptp = xlim[0] - xlim[1]
pad = (ptp / np.cos(np.mean(ylim)) - ptp) / 2
xlim = [xlim[0] + pad, xlim[1] - pad]

plt.setp(ax, xlim=xlim)
plt.setp(
    ax,
    xlabel="RA (deg)",
    ylim=ylim,
    ylabel="Dec (deg)",
    aspect=1,
    adjustable="datalim",
)
# plt.setp(ax, xlim=ax.get_xlim()[::-1])
ax.minorticks_on()
ax.grid()
plt.tight_layout()

if args.o:
    plt.savefig(args.o, dpi=200)

plt.show()
