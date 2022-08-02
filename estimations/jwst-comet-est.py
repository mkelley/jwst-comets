#!/usr/bin/env python3
# Licensed with the MIT License, see LICENSE.txt
# Author: Michael S. P. Kelley
import sys
import argparse
from collections import OrderedDict
import logging

import numpy as np
import astropy.units as u
from astropy.table import Table
from sbpy.activity import Afrho, Efrho
from sbpy.activity import RectangularAperture


def aperture(a):
    if a == "ifu":
        aper = u.Quantity(0.2, u.arcsec)
    elif a == "fs200":
        aper = RectangularAperture(u.Quantity((0.2, 1.0), u.arcsec))
    elif a == "fs400":
        aper = RectangularAperture(u.Quantity((0.4, 1.2), u.arcsec))
    elif a == "fs1600":
        aper = RectangularAperture(u.Quantity((1.6, 1.6), u.arcsec))
    elif a == "lrs":
        aper = RectangularAperture(u.Quantity((0.5, 1.6), u.arcsec))
    elif a == "mrs":
        aper = u.Quantity(1.3, u.arcsec)
    elif a == "nircam":
        aper = u.Quantity(0.13, u.arcsec)
    else:
        aper = u.Quantity(float(a), u.arcsec)
    # else:
    #    raise argparse.ArgumentError('Invalid --aper: {}'.format(str(a)))

    return aper


def constant_spectral_resolution(start, stop, R):
    # from mskpy
    d = 1 + 1 / R
    n = int(np.ceil(np.log(stop / start) / np.log(d)))
    return start * d ** np.arange(n)


def estimate(args):
    """args - command-line arguments"""
    # deduce phase angle, assume JWST at 1.0 au
    cosPhase = (args.rh**2 + args.delta**2 - 1.0**2) / (
        2 * args.rh * args.delta
    )
    phase = np.degrees(np.arccos(cosPhase))
    logging.info("Phase angle = {:.1f}".format(phase))
    eph = dict(
        rh=u.Quantity(args.rh, u.au),
        delta=u.Quantity(args.delta, u.au),
        phase=u.Quantity(phase, u.deg),
    )

    meta = OrderedDict()
    meta["cmd"] = " ".join(sys.argv)
    meta["rh"] = str(eph["rh"])
    meta["delta"] = str(eph["delta"])
    meta["phase"] = str(eph["phase"])
    meta["aper"] = str(args.aper)

    if args.afrho is not None:
        est = dust_estimate(eph, args, meta)
    else:
        raise NotImplementedError("gas estimates are not yet implemented")
        # est = gas_estimate(eph, args, meta)

    return est


def dust_estimate(eph, args, meta):
    wave = constant_spectral_resolution(0.5, 30, args.R)
    wave = u.Quantity(wave, "um")
    if args.m:
        mH = args.m - 5 * np.log10(eph["delta"].to("au").value)
        afrho = Afrho(10 ** (-0.17 * mH + 3.77), "cm")
        afrho *= 10 ** (0.69 * (args.dusty - args.gassy))
    else:
        afrho = Afrho(args.afrho, "cm")

    efrho = Efrho(args.ef2af * afrho)

    if isinstance(args.aper, RectangularAperture):
        area = args.aper.shape[0] * args.aper.shape[1]
    else:
        area = np.pi * args.aper**2

    if args.unit.is_equivalent("mJy/arcsec2") or args.unit.is_equivalent(
        "W/(m2 um sr)"
    ):
        # user requested surface brightness units
        surface_brightness = True
        unit = args.unit * area
    else:
        surface_brightness = False
        unit = args.unit

    fsca = afrho.to_fluxd(wave, args.aper, eph, phasecor=True, unit=unit)
    fth = efrho.to_fluxd(wave, args.aper, eph, unit=unit, Tscale=args.Tscale)
    if surface_brightness:
        fsca = (fsca / area).to(args.unit)
        fth = (fth / area).to(args.unit)

    est = Table(
        data=(wave, fsca + fth, fsca, fth),
        names=("wave", "total", "F_sca", "F_th"),
    )

    est.meta = meta
    est.meta["phase function"] = float(
        Afrho(1, "cm").to_phase(eph["phase"], 0 * u.deg).value
    )
    est.meta["ef2af"] = args.ef2af
    est.meta["Tscale"] = args.Tscale
    est.meta["Afrho"] = str(afrho)
    est.meta["efrho"] = str(efrho)

    return est


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate model comet spectra for JWST ETC.",
        epilog="For gas output, -m uses the Jorda et al. (2008, ACM, 8046) correlation Q=10**(30.675 - 0.2453 * mH), where mH is heliocentric magnitude.  This assumes --h2o=1.0.  --gassy and --dusty can be repeatedly used to scale the results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("rh", type=float, help="heliocentric distance (au)")
    parser.add_argument(
        "delta", type=float, help="observer-target distance (au)"
    )

    parser.add_argument(
        "--aper",
        type=aperture,
        default=aperture(0.2),
        help="aperture radius (arcsec), or one of: ifu, fs200, fs400, fs1600, lrs, mrs, nircam",
    )

    parser.add_argument("-R", default=10, type=int, help="spectral resolution")

    parser.add_argument(
        "--unit",
        type=u.Unit,
        default="mJy/arcsec2",
        help="output spectral flux density or surface brightness unit",
    )

    parser.add_argument(
        "-m",
        action="store_true",
        help="indicates that --afrho or -Q is a total visual mangitude that should be converted before using",
    )

    parser.add_argument(
        "--dusty",
        "-d",
        action="count",
        default=0,
        help="indicates that the conversion from -m should assume a dusty coma, repeated options increase dust/gas: Q/1.6|2.4|3.7, Afrho*4.8|23|110",
    )

    parser.add_argument(
        "--gassy",
        "-g",
        action="count",
        default=0,
        help="indicates that the conversion from -m should assume a gassy coma, repeated options increase gas/dust: Q*1.6|2.4|3.7, Afrho/4.8|23|110",
    )

    dust_or_gas = parser.add_mutually_exclusive_group(required=True)
    dust_or_gas.add_argument(
        "--afrho",
        type=float,
        help="dust mode based on this Afρ at 0 deg phase (cm) or total visual magnitude (EXPERIMENTAL, requires -m)",
    )
    dust_or_gas.add_argument(
        "-Q",
        type=float,
        help="gas mode based on this scale factor (molecules/s) or total visual magnitude (requires -m)",
    )

    dust = parser.add_argument_group(title="dust options")
    dust.add_argument(
        "--ef2af", type=float, default=3.5, help="the ratio εfρ/Afρ"
    )
    dust.add_argument(
        "--Tscale",
        type=float,
        default=1.1,
        help="LTE blackbody temperature scale factor",
    )

    gas = parser.add_argument_group(title="gas options")
    gas.add_argument(
        "--h2o", type=float, default=1.0, help="relative H2O production rate"
    )
    gas.add_argument(
        "--co2", type=float, default=0.15, help="relative CO2 production rate"
    )
    gas.add_argument(
        "--co", type=float, default=0.05, help="relative CO production rate"
    )

    parser.add_argument("-o", help="save to this file name")

    args = parser.parse_args()
    est = estimate(args)
    outf = sys.stdout if args.o is None else args.o
    est.write(outf, format="ascii.ecsv")
