import numpy as np
import astropy.units as u
from sbpy.activity import Afrho

# Goal is to reproduce flux density within a 62 mas radius aperture (~2
# NIRCam SWC pixels).
afrho = Afrho(500 * u.cm)
filt = 'NIRCAM F182M'
rap = 62 * u.mas
eph = {  # ignoring phase effects
    'rh': 1.52 * u.au,
    'delta': 2.34 * u.au
}
fluxd0 = afrho.to_fluxd(filt, rap, eph, unit=u.mJy)
print('Goal flux density: {:.4f}'.format(fluxd0))
# 0.0852 mJy

# scale is 0.1 mas / pix
y, x = np.mgrid[-1000:1000, -1000:1000] / 10
rho = np.hypot(y, x)

rho[rho <= 1] = 1
coma = 1 / rho

# ETC min core size is 1 mas.  What is the mean surface brightness of a coma
# within 1 mas?
fluxd1 = afrho.to_fluxd(filt, 1 * u.mas, eph, unit=u.mJy)
sb1 = (fluxd1 / (1 * u.mas)**2 / np.pi)

# Scale the coma to this surface brightness
coma = coma * sb1

# Resample image to 1 mas / pix, compare with goal flux density.
coma_rescaled = (
    np.mean(np.mean(coma.reshape(200, 10, 2000), 1).reshape(200, 200, 10), 2))
rho_rescaled = (
    np.mean(np.mean(rho.reshape(200, 10, 2000), 1).reshape(200, 200, 10), 2))

aperture = rho_rescaled <= rap.to_value('mas')
fluxd = coma_rescaled[aperture].sum() * (1 * u.mas)**2
print('{} aperture from image (and ratio): {:.4f} ({:.4f})'.format(
    rap, fluxd, fluxd / fluxd0))
# 0.1690 mJy (1.9840)

# The result is too bright by a factor of 2.  This is because the ETC profile
# scales from the disk's mean surface brightness, but a true coma scales from
# the center of disk.  It can be shown that the mean surface brightness of a
# coma within an aperture is the surface brightness at 1/2 the aperture radius.
# Repeat the exercise, but with the disk at 1/2 the brightness (using 0.5 mas
# radius in the Afrho call).

# NOTE THAT THE SURFACE BRIGHTNESS IS COMPUTED WITH 1 mas!
fluxd = afrho.to_fluxd(filt, 0.5 * u.mas, eph, unit=u.mJy)
coma2 = 1 / rho * fluxd / (1 * u.mas)**2 / np.pi
coma_rescaled2 = (
    np.mean(np.mean(coma2.reshape(200, 10, 2000), 1).reshape(200, 200, 10), 2))

fluxd = coma_rescaled2[aperture].sum() * (1 * u.mas)**2
print('{} aperture from revised image (and ratio): {:.4f} ({:.4f})'.format(
    rap, fluxd, fluxd / fluxd0))
# 0.0845 mJy (0.9920)

# The result is low by 0.8%, an acceptable tolerance for ETC purposes.
