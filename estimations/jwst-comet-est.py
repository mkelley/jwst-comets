import argparse
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

def estimate(args):
    """args - command-line arguments"""
    from sbpy.activity import Afrho, Efrho

    # deduce phase angle, assume JWST at 1.0 au
    cosPhase = ((args.rh**2 + args.delta**2 - 1.0**2)
                / (2 * args.rh * args.delta))
    eph = dict(
        rh=u.Quantity(args.rh, u.au),
        delta=u.Quantity(args.delta, u.au),
        phase=u.Quantity(cosPhase, u.rad)
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate model comet spectra for JWST ETC.',
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
        help='aperture radius (arcsec), or one of: ifu, fs200, fs400, fs1600, lrs, mrs, nircam')

    dust_or_gas = parser.add_mutually_exclusive_group(required=True)
    dust_or_gas.add_argument(
        '--afrho',
        type=float,
        help='dust mode based on this Afρ value (cm)')
    dust_or_gas.add_argument(
        '-Q',
        type=float,
        help='gas mode based on this scale factor (molecules/s)')

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
    estimate(args)
