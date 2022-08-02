# APT "fixed" moving targets

Author: Michael S. P. Kelley

## Background

As of May 2021, APT does not visualize moving targets. Instead, fixed targets may be used as a work around. The best fixed targets use the ephemeris of the moving target in question. To that end, the script `apt-fixed-moving-target.py` generates a set of fixed targets centered on the object ephemeris.  You specify the date range and step size.  Only times when the target is observable by JWST are returned.


## Requirements
* Python 3
* astropy
* astroquery

## Usage
This is a command-line Python script.  A brief help may be displayed:

```
$ python3 apt-fixed-moving-target.py -h
usage: apt-fixed-moving-target [-h] [--type {comet,designation,smallbody,majorbody}] [--step STEP] [-o O] [-f] [--xml] [--nircam] [--no-cache] target start_date stop_date

Generate fixed targets for JWST APT based on moving target ephemerides.

positional arguments:
  target                moving target, e.g., 1P, 24, P/2003 S2
  start_date            UTC start date, YYYY-MM-DD
  stop_date             UTC stop date, YYYY-MM-DD

optional arguments:
  -h, --help            show this help message and exit
  --type {comet,designation,smallbody,majorbody}
                        target type
  --step STEP           Step size, e.g., 5d
  -o O                  output to this file name
  -f                    force overwrite output file
  --xml                 format output as APT XML
  --nircam              observe with NIRCam (for XML output)
  --no-cache            do not use cached ephemeris

Target names must be resolvable by JPL Horizons. Specifying --type=comet will use the "closest apparition" and "no fragment" search flags.
```

There are two output modes: a plain text target list, and a JWST XML APT file.
* The plain text target list can be imported into APT and each fixed target
assigned to whatever observations you want.
* The XML file can optionally include some dummy NIRCam observations,
one for each ephemeris epoch.  The observations have position angle
constraints, so that when they are viewed with APT, they already have
the approximate (to the best of my knowledge) orientation for that
ephemeris epoch.

## Examples
### Plain text table output
A plain text table of targets can be generated.  This list helps identify the time periods an object is observable, and the position in the sky when it is observable.  It does not help understand the orientation of the instruments at those times.

Generate fixed targets for comet 238P, output as space-separated-value table
```
$ python3 apt-fixed-moving-target.py 116P 2022-05-01 2022-07-01 --type=comet -o comet.txt
$ cat comet.txt
# Name RA DEC "RA Uncertainty" "Dec Uncertainty" Comments
116P.Wild.4.at.2022-05-01T00 144.61541 16.73492 2.456 0.968 "2022-05-01 00:00:00.000 UTC, rh: 2.262 au, delta: 1.83 au, phase: 26.18 deg, target->Sun PA: 290 deg, velocity PA: 291 deg"
116P.Wild.4.at.2022-05-06T00 145.78329 16.14966 2.397 0.965 "2022-05-06 00:00:00.000 UTC, rh: 2.254 au, delta: 1.88 au, phase: 26.61 deg, target->Sun PA: 290 deg, velocity PA: 292 deg"
116P.Wild.4.at.2022-05-11T00 147.0766 15.5206 2.342 0.963 "2022-05-11 00:00:00.000 UTC, rh: 2.247 au, delta: 1.93 au, phase: 26.91 deg, target->Sun PA: 290 deg, velocity PA: 292 deg"
116P.Wild.4.at.2022-05-16T00 148.48423 14.84954 2.291 0.963 "2022-05-16 00:00:00.000 UTC, rh: 2.240 au, delta: 1.98 au, phase: 27.10 deg, target->Sun PA: 290 deg, velocity PA: 293 deg"
116P.Wild.4.at.2022-05-21T00 149.99567 14.13821 2.245 0.965 "2022-05-21 00:00:00.000 UTC, rh: 2.233 au, delta: 2.03 au, phase: 27.19 deg, target->Sun PA: 290 deg, velocity PA: 293 deg"
```

Same, but for an asteroid and a planet:
```
$ python3 apt-fixed-moving-target.py 24 2022-05-01 2023-05-01 -o asteroid.txt
$ python3 apt-fixed-moving-target.py 599 2022-05-01 2023-05-01 --type=majorbody -o planet.txt
```

See JDOX for instructions on how to [read these target lists](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview/apt-workflow-articles/apt-bulk-target-ingest#APTBulkTargetIngest-AccessingtheImporter) into APT.

### XML output and generated observations
Use `--xml` to generate an APT file in XML format, which can be imported into APT (File->Import->JWST XML file).

Templated observations of your target may also be generated for NIRCam using `--nircam`.  The observations will have position angle and timing constraints based on the target-Sun angle and ephemeris epoch.  With the position angle constraint, the Aladin viewer will automatically give NIRCam the approximate orientation.  These constraints should not be used for designing real observations.

---
**Warning**

The orientations are based on the target-Sun angle, which, to the best of my knowledge, is correct.  However, please verify the orientation yourself and let me know of any successes or failures.

---

```
python3 apt-fixed-moving-target.py 116P 2022-05-01 2022-07-01 --type=comet --xml --nircam -o comet.xml
```
![APT, example NIRCam observation](apt-comet-example-vis.png)
![APT, multiple example NIRCam observations](apt-comet-multi-example-vis.png)