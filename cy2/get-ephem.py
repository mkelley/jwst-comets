import os
import numpy as np
from astropy.io import ascii
from astroquery.mpc import MPC


targets = ascii.read("ztf.txt")
mpc = MPC()

errored = []
for target in targets:
    if target["designation"] is np.ma.masked:
        continue

    outf = "{}.csv".format(
        target["designation"].replace(" ", "").replace("/", "").lower()
    )
    if os.path.exists(outf):
        continue
    try:
        eph = mpc.get_ephemeris(
            target["designation"],
            start="2023-07-01",
            step="5d",
            number=300,
            unc_links=False,
            cache=True,
        )
    except:
        errored.append(target["designation"])
        continue

    i = (eph["Elongation"] > 85) * (eph["Elongation"] < 135)
    eph[i].write(outf)

print("errored:\n", "\n".join(errored))
