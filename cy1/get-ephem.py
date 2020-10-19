import os
from astropy.io import ascii
from astroquery.mpc import MPC


targets = ascii.read('ztf.txt')
mpc = MPC()

errored = []
for target in targets:
    outf = '{}.csv'.format(
        target['desg'].replace(' ', '').replace('/', '').lower()
    )
    if os.path.exists(outf):
        continue
    try:
        eph = mpc.get_ephemeris(target['desg'], start='2022-05-01', step='5d',
                                number=150, unc_links=False, cache=True)
    except:
        errored.append(target['desg'])
        continue

    i = (eph['Elongation'] > 85) * (eph['Elongation'] < 135)
    eph[i].write(outf)

print('errored:\n', '\n'.join(errored))
