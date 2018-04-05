import sys
import argparse
import numpy as np
import astropy.units as u

def aperture(a):
    from sbpy.activity import RectangularAperture
    
    if a == 'ifu':
        aper = u.Quantity(0.2, u.arcsec)
    elif a == 'fs200':
        aper = RectangularAperture(u.Quantity((0.2, 1.0), u.arcsec))
    elif a == 'fs400':
        aper = RectangularAperture(u.Quantity((0.4, 1.2), u.arcsec))
    elif a == 'fs1600':
        aper = RectangularAperture(u.Quantity((1.6, 1.6), u.arcsec))
    elif a == 'lrs':
        aper = RectangularAperture(u.Quantity((0.5, 1.6), u.arcsec))
    elif a == 'mrs':
        aper = u.Quantity(1.3, u.arcsec)
    elif a == 'nircam':
        aper = u.Quantity(0.13, u.arcsec)
    elif isinstance(a, (float, int)):
        aper = u.Quantity(a, u.arcsec)
    else:
        raise argparse.ArgumentError('Invalid --aper: {}'.format(str(aper)))

    return aper

def constant_spectral_resolution(start, stop, R):
    # from mskpy
    d = 1 + 1 / R
    n = int(np.ceil(np.log(stop / start) / np.log(d)))
    return start * d**np.arange(n)

def estimate(args):
    """args - command-line arguments"""
    from collections import OrderedDict
    
    # deduce phase angle, assume JWST at 1.0 au
    cosPhase = ((args.rh**2 + args.delta**2 - 1.0**2)
                / (2 * args.rh * args.delta))
    eph = dict(
        rh=u.Quantity(args.rh, u.au),
        delta=u.Quantity(args.delta, u.au),
        phase=u.Quantity(np.degrees(np.arccos(cosPhase)), u.deg)
    )

    meta = OrderedDict()
    meta['cmd'] = ' '.join(sys.argv)
    meta['rh'] = str(eph['rh'])
    meta['delta'] = str(eph['delta'])
    meta['phase'] = str(eph['phase'])
    meta['aper'] = str(args.aper)

    if args.afrho is not None:
        est = dust_estimate(eph, args, meta)
    else:
        raise NotImplemented('gas estimates are not yet implemented')
        est = gas_estimate(eph, args, meta)

    return est

def dust_estimate(eph, args, meta):
    from astropy.table import Table
    from sbpy.activity import Afrho, Efrho

    wave = constant_spectral_resolution(0.5, 30, args.R)
    if args.m:
        mH = args.m - 5 * np.log10(eph['delta'])
        afrho = Afrho(10**(-0.21 * mH + 5.69), 'cm')

        if args.gassy:
            afrho /= 2.9
        elif args.dusty:
            afrho *= 2.9
    else:
        afrho = Afrho(args.afrho, 'cm')

    efrho = Efrho(args.ef2af * afrho)

    fsca = afrho.fluxd(wave, args.aper, eph, phasecor=True, unit='mJy')
    fth = efrho.fluxd(wave, args.aper, eph, unit='mJy', Tscale=args.Tscale)
    est = Table(data=(wave, fsca + fth, fsca, fth),
                names=('wave', 'total', 'F_sca', 'F_th'))

    est.meta = meta
    est.meta['ef2af'] = args.ef2af
    est.meta['Tscale'] = args.Tscale
    est.meta['Afrho'] = str(afrho)
    est.meta['efrho'] = str(efrho)

    return est

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate model comet spectra for JWST ETC.',
        epilog='For gas output, -m uses the Jorda et al. (2008, ACM, 8046) correlation Q=10**(30.675 - 0.2453 * mH), where mH is heliocentric magnitude.  This assumes --h2o=1.0.  --gassy scales Q up a factor of 1.6, --dusty scales Q down by a factor of 1.6.  For dust output, -m uses MSK\'s prototype correlation 10**(-0.21 * mH + 5.69), scaled by 2.9 for --dusty and --gassy.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'rh',
        type=float,
        help='heliocentric distance (au)')
    parser.add_argument(
        'delta',
        type=float,
        help='observer-target distance (au)')

    parser.add_argument(
        '--aper',
        type=aperture,
        default=aperture(0.2),
        help='aperture radius (arcsec), or one of: ifu, fs200, fs400, fs1600, lrs, mrs, nircam')

    parser.add_argument(
        '-R',
        default=10,
        help='spectral resolution'
    )

    parser.add_argument(
        '-m',
        action='store_true',
        help='indicates that --afrho or -Q is a total visual mangitude that should be converted before using')

    parser.add_argument(
        '--dusty',
        action='store_true',
        help='indicates that the conversion from -m should assume a dusty coma')

    parser.add_argument(
        '--gassy',
        action='store_true',
        help='indicates that the conversion from -m should assume a gassy coma')

    dust_or_gas = parser.add_mutually_exclusive_group(required=True)
    dust_or_gas.add_argument(
        '--afrho',
        type=float,
        help='dust mode based on this Afρ at 0 deg phase (cm) or total visual magnitude (requires -m)')
    dust_or_gas.add_argument(
        '-Q',
        type=float,
        help='gas mode based on this scale factor (molecules/s) or total visual magnitude (requires -m)')

    dust = parser.add_argument_group(title='dust options')
    dust.add_argument(
        '--ef2af',
        type=float,
        default=3.5,
        help='the ratio εfρ/Afρ')
    dust.add_argument(
        '--Tscale',
        type=float,
        default=1.1,
        help='LTE blackbody temperature scale factor')

    gas = parser.add_argument_group(title='gas options')
    gas.add_argument(
        '--h2o',
        type=float,
        default=1.0,
        help='relative H2O production rate')
    gas.add_argument(
        '--co2',
        type=float,
        default=0.15,
        help='relative CO2 production rate')
    gas.add_argument(
        '--co',
        type=float,
        default=0.05,
        help='relative CO production rate')

    args = parser.parse_args()
    est = estimate(args)
    est.write(sys.stdout, format='ascii.ecsv')
