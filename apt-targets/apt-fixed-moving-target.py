#!/usr/bin/env python3
# Licensed with the MIT License, see LICENSE.txt
# Author: Michael S. P. Kelley
import os
import sys
import argparse

try:
    from lxml import etree
except ImportError:
    etree = None

import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import Angle
from astroquery.jplhorizons import Horizons

parser = argparse.ArgumentParser(
    "apt-fixed-moving-target",
    description=(
        "Generate fixed targets for JWST APT based on moving target"
        " ephemerides."
    ),
    epilog="""Target names must be resolvable by JPL Horizons.
Specifying --type=comet will use the "closest apparition" and "no fragment" search flags.""",
)
parser.add_argument("target", help="moving target, e.g., 1P, 24, P/2003 S2")
parser.add_argument("start_date", type=Time, help="UTC start date, YYYY-MM-DD")
parser.add_argument("stop_date", type=Time, help="UTC stop date, YYYY-MM-DD")
parser.add_argument(
    "--type",
    choices=["comet", "designation", "smallbody", "majorbody"],
    default="smallbody",
    help="target type",
)
parser.add_argument("--step", default="5d", help="Step size, e.g., 5d")
parser.add_argument("-o", default=sys.stdout, help="output to this file name")
parser.add_argument(
    "-f", action="store_true", help="force overwrite output file"
)
parser.add_argument(
    "--xml", action="store_true", help="format output as APT XML"
)
parser.add_argument(
    "--nircam", action="store_true", help="observe with NIRCam (for XML output)"
)
parser.add_argument(
    "--nirspec-ifu",
    action="store_true",
    help="observe with NIRSpec IFU (for XML output)",
)
parser.add_argument(
    "--no-cache", action="store_true", help="do not use cached ephemeris"
)
args = parser.parse_args()

# parameter checks
if isinstance(args.o, str) and os.path.exists(args.o) and not args.f:
    print("File exists, use -f to force overwrite.")
    sys.exit(1)

if args.xml and etree is None:
    print("XML output requires lxml.")
    sys.exit(1)


def get_ephemeris(args):
    epochs = dict(
        start=args.start_date.iso, stop=args.stop_date.iso, step=args.step
    )
    id_type = args.type

    # ephemeris options
    opts = dict(solar_elongation=(85, 135), cache=(not args.no_cache))

    # comets need special care: use closest orbital elements by date, do not
    # search for related fragments
    if args.type == "comet":
        id_type = "designation"
        opts.update(closest_apparition=True, no_fragments=True)

    horizons = Horizons(
        args.target, location="@jwst", epochs=epochs, id_type=id_type
    )
    eph = horizons.ephemerides(**opts)
    eph["date"] = Time(eph["datetime_jd"], format="jd")

    if len(eph) == 0:
        print(
            """No epochs observable by JWST.
    Try a new date range and/or finer time steps."""
        )
        sys.exit(1)

    return eph


def format_comments(row):
    "Target comments."
    return (
        f"""{row['date'].iso} UTC, """
        f"""rh: {row['r']:.3f} au, """
        f"""delta: {row['delta']:.3g} au, """
        f"""phase: {row['alpha']:.2f} deg, """
        f"""target->Sun PA: {(Angle(row['sunTargetPA'], 'deg') - 180 * u.deg).wrap_at(360 * u.deg):.0f}, """
        f"""velocity PA: {(Angle(row['velocityPA'], 'deg')).wrap_at(360 * u.deg):.0f}"""
    )


def format_between(row):
    """Observation between dates contraint.

    For example:
      After="{17-SEP-2022:00:00:00}" Before="19-SEP-2022:00:00:00"

    """

    after = (
        (row["date"] - 12 * u.hr).datetime.strftime("%d-%b-%Y:%H:00:00").upper()
    )
    before = (
        (row["date"] + 12 * u.hr).datetime.strftime("%d-%b-%Y:%H:00:00").upper()
    )
    return after, before


def eph_to_table(eph):
    """Format ephemeris as an APT target table.

    https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview/apt-workflow-articles/apt-bulk-target-ingest

    Columns: Name, RA, DEC, RA Uncertainty, Dec Uncertainty, RA PM, RA PM units, DEC PM, DEC PM units, Epoch, Annual Parallax, Comments
    """

    tab = eph["RA", "DEC", "RA_3sigma", "DEC_3sigma"]
    tab["RA_3sigma"].name = "RA Uncertainty"
    tab["DEC_3sigma"].name = "Dec Uncertainty"
    tab["Comments"] = [format_comments(row) for row in eph]
    target = "".join([c if c.isalnum() else "." for c in eph["targetname"][0]])
    tab["Name"] = [f"""{target}.at.{row['date'].isot[:13]}""" for row in eph]

    # reorder columns
    tab = tab[
        "Name", "RA", "DEC", "RA Uncertainty", "Dec Uncertainty", "Comments"
    ]
    return tab


def eph_to_xml(eph):
    """Format ephemeris as APT XML file."""
    xml_parser = etree.XMLParser(remove_blank_text=True)

    root = etree.XML(
        f"""<JwstProposal schemaVersion="51" APTVersion="Version 2020.5  " PRDVersion="PRDOPSSOC-032" xmlns="http://www.stsci.edu/JWST/APT" xmlns:mmrscge="http://www.stsci.edu/JWST/APT/Template/MiriMRSCrossGratingEngineering" xmlns:nsmsasd="http://www.stsci.edu/JWST/APT/Template/NirspecMSAShortDetect" xmlns:nid="http://www.stsci.edu/JWST/APT/Template/NirissDark" xmlns:ncipri="http://www.stsci.edu/JWST/APT/Template/NircamIprImaging" xmlns:nif="http://www.stsci.edu/JWST/APT/Template/NirissFocus" xmlns:wfscfp="http://www.stsci.edu/JWST/APT/Template/WfscFinePhasing" xmlns:nii="http://www.stsci.edu/JWST/APT/Template/NirissImaging" xmlns:mmimf="http://www.stsci.edu/JWST/APT/Template/MiriMimf" xmlns:mlrs="http://www.stsci.edu/JWST/APT/Template/MiriLRS" xmlns:mc="http://www.stsci.edu/JWST/APT/Template/MiriCoron" xmlns:md="http://www.stsci.edu/JWST/APT/Template/MiriDark" xmlns:nsfgwt="http://www.stsci.edu/JWST/APT/Template/NirspecFilterGratingWheelTest" xmlns:wfscga="http://www.stsci.edu/JWST/APT/Template/WfscGlobalAlignment" xmlns:nsil="http://www.stsci.edu/JWST/APT/Template/NirspecInternalLamp" xmlns:idfu="http://www.stsci.edu/JWST/APT/Template/IsimDictionaryFileUpdate" xmlns:mi="http://www.stsci.edu/JWST/APT/Template/MiriImaging" xmlns:nsimg="http://www.stsci.edu/JWST/APT/Template/NirspecImaging" xmlns:niami="http://www.stsci.edu/JWST/APT/Template/NirissAmi" xmlns:niwfss="http://www.stsci.edu/JWST/APT/Template/NirissWfss" xmlns:ncif="http://www.stsci.edu/JWST/APT/Template/NircamInternalFlat" xmlns:nsfss="http://www.stsci.edu/JWST/APT/Template/NirspecFixedSlitSpectroscopy" xmlns:fgsif="http://www.stsci.edu/JWST/APT/Template/FgsInternalFlat" xmlns:ncef="http://www.stsci.edu/JWST/APT/Template/NircamExternalFlat" xmlns:ncei="http://www.stsci.edu/JWST/APT/Template/NircamEngineeringImaging" xmlns:niec="http://www.stsci.edu/JWST/APT/Template/NirissExternalCalibration" xmlns:nsmsam="http://www.stsci.edu/JWST/APT/Template/NirspecMSAMasking" xmlns:niif="http://www.stsci.edu/JWST/APT/Template/NirissInternalFlat" xmlns:ncpili="http://www.stsci.edu/JWST/APT/Template/NircamPilImaging" xmlns:nsmsaa="http://www.stsci.edu/JWST/APT/Template/NirspecMSAAnneal" xmlns:mcpc="http://www.stsci.edu/JWST/APT/Template/MiriCpc" xmlns:fgsec="http://www.stsci.edu/JWST/APT/Template/FgsExternalCalibration" xmlns:nsd="http://www.stsci.edu/JWST/APT/Template/NirspecDark" xmlns:nsf="http://www.stsci.edu/JWST/APT/Template/NirspecFocus" xmlns:ns="http://www.stsci.edu/JWST/APT/Instrument/Nirspec" xmlns:nisoss="http://www.stsci.edu/JWST/APT/Template/NirissSoss" xmlns:nsmimf="http://www.stsci.edu/JWST/APT/Template/NirspecMimf" xmlns:ncts="http://www.stsci.edu/JWST/APT/Template/NircamTimeSeries" xmlns:ncd="http://www.stsci.edu/JWST/APT/Template/NircamDark" xmlns:mef="http://www.stsci.edu/JWST/APT/Template/MiriExternalFlat" xmlns:ncc="http://www.stsci.edu/JWST/APT/Template/NircamCoron" xmlns:ncf="http://www.stsci.edu/JWST/APT/Template/NircamFocus" xmlns:mmrs="http://www.stsci.edu/JWST/APT/Template/MiriMRS" xmlns:nci="http://www.stsci.edu/JWST/APT/Template/NircamImaging" xmlns:sk="http://www.stsci.edu/JWST/APT/Template/StationKeeping" xmlns:wfscc="http://www.stsci.edu/JWST/APT/Template/WfscCommissioning" xmlns:iat="http://www.stsci.edu/JWST/APT/Template/IsimAsicTuning" xmlns:sr="http://www.stsci.edu/JWST/APT/Template/SafeModeRecovery" xmlns:rtc="http://www.stsci.edu/JWST/APT/Template/RealtimeCommanding" xmlns:nsfr="http://www.stsci.edu/JWST/APT/Template/NirspecFocusReference" xmlns:ncw="http://www.stsci.edu/JWST/APT/Template/NircamWheelExercise" xmlns:nsbots="http://www.stsci.edu/JWST/APT/Template/NirspecBrightObjectTimeSeries" xmlns:mann="http://www.stsci.edu/JWST/APT/Template/MiriAnneal" xmlns:nsmos="http://www.stsci.edu/JWST/APT/Template/NirspecMOS" xmlns:wfscco="http://www.stsci.edu/JWST/APT/Template/WfscControlOnly" xmlns:ncgts="http://www.stsci.edu/JWST/APT/Template/NircamGrismTimeSeries" xmlns:wfsccp="http://www.stsci.edu/JWST/APT/Template/WfscCoarsePhasing" xmlns:msa="http://www.stsci.edu/JWST/APT/Template/NirspecMSA" xmlns:nsifus="http://www.stsci.edu/JWST/APT/Template/NirspecIFUSpectroscopy" xmlns:ncwfss="http://www.stsci.edu/JWST/APT/Template/NircamWfss" xmlns:nifs="http://www.stsci.edu/JWST/APT/Template/NirissFlatSuite" xmlns:fgsf="http://www.stsci.edu/JWST/APT/Template/FgsFocus" xmlns:po="http://www.stsci.edu/JWST/APT/Template/PointingOnly">
    <ProposalInformation>
        <ProposalPhase>Draft</ProposalPhase>
        <StsciEditNumber>0</StsciEditNumber>
        <ProposalCategory>GO</ProposalCategory>
        <ProposalSize>SMALL</ProposalSize>
        <ProprietaryPeriod>C[12 Months]</ProprietaryPeriod>
        <Cycle>1</Cycle>
        <ParallelNextCycle>false</ParallelNextCycle>
        <ParallelThirdCycle>false</ParallelThirdCycle>
        <NoSimilarSubmissions>false</NoSimilarSubmissions>
        <PrincipalInvestigator>
            <InvestigatorAddress>
                <ESAMember>false</ESAMember>
                <CSAMember>false</CSAMember>
                <Retired>false</Retired>
            </InvestigatorAddress>
            <Contact>true</Contact>
            <GTOTimeSource></GTOTimeSource>
        </PrincipalInvestigator>
        <CoInvestigator>
            <CoPI>false</CoPI>
            <InvestigatorAddress>
                <ESAMember>false</ESAMember>
                <CSAMember>false</CSAMember>
                <Retired>false</Retired>
            </InvestigatorAddress>
            <AdminUSPI>false</AdminUSPI>
            <Contact>false</Contact>
        </CoInvestigator>
    </ProposalInformation>
    <Targets></Targets>
    <PureParallelSlots/>
    <DataRequests/>
    <LinkingRequirements/>
</JwstProposal>
    """,
        xml_parser,
    )

    targets = root.find("Targets", {None: "http://www.stsci.edu/JWST/APT"})
    data_requests = root.find(
        "DataRequests", {None: "http://www.stsci.edu/JWST/APT"}
    )

    if args.nircam:
        nircam = etree.XML(
            """
<ObservationGroup>
    <Label>NIRCam</Label>
</ObservationGroup>
""",
            xml_parser,
        )

    if args.nirspec_ifu:
        nirspec_ifu = etree.XML(
            """
<ObservationGroup>
    <Label>NIRSpec IFU</Label>
</ObservationGroup>
""",
            xml_parser,
        )

    target = "".join([c if c.isalnum() else "." for c in eph["targetname"][0]])
    for i in range(len(eph)):
        target_id = f"{target}.at.{eph['date'][i].isot[:13]}"
        comments = format_comments(eph[i])
        targets.append(
            etree.XML(
                f"""
    <Target xsi:type="FixedTargetType" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
        <Number>{i + 1}</Number>
        <TargetName>{target_id}</TargetName>
        <TargetID>{target_id}</TargetID>
        <Comments>{comments}</Comments>
        <RAProperMotion></RAProperMotion>
        <DecProperMotion></DecProperMotion>
        <RAProperMotionUnits></RAProperMotionUnits>
        <DecProperMotionUnits></DecProperMotionUnits>
        <Extended>Unknown</Extended>
        <EquatorialCoordinates Value="{Angle(eph['RA'][i], 'deg').to_string(unit='hr', sep=' ')} {Angle(eph['DEC'][i], 'deg').to_string(unit='deg', sep=' ', alwayssign=True)}">
            <RAUncertainty>
                <Value>{eph['RA_3sigma'][i]}</Value>
                <Units>Arcsec</Units>
            </RAUncertainty>
            <DecUncertainty>
                <Value>{eph['DEC_3sigma'][i]}</Value>
                <Units>Arcsec</Units>
            </DecUncertainty>
        </EquatorialCoordinates>
        <BackgroundTargetReq>false</BackgroundTargetReq>
        <TargetConfirmationRun>false</TargetConfirmationRun>
    </Target>
""",
                xml_parser,
            )
        )

        if args.nircam:
            after, before = format_between(eph[i])
            nircam.append(
                etree.XML(
                    f"""
<Observation AutoTarget="false">
    <Number>{i + 1}</Number>
    <TargetID>{target_id}</TargetID>
    <Instrument>NIRCAM</Instrument>
    <Template xmlns:nci="http://www.stsci.edu/JWST/APT/Template/NircamImaging">
        <nci:NircamImaging>
            <nci:Module>ALL</nci:Module>
            <nci:Subarray>FULL</nci:Subarray>
            <nci:PrimaryDitherType>INTRAMODULEX</nci:PrimaryDitherType>
            <nci:PrimaryDithers>5</nci:PrimaryDithers>
            <nci:SubpixelDitherType>STANDARD</nci:SubpixelDitherType>
            <nci:SubpixelPositions>1</nci:SubpixelPositions>
            <nci:Filters>
                <nci:FilterConfig>
                    <nci:ShortFilter>F070W</nci:ShortFilter>
                    <nci:LongFilter>F277W</nci:LongFilter>
                    <nci:ReadoutPattern>BRIGHT1</nci:ReadoutPattern>
                    <nci:Groups>10</nci:Groups>
                    <nci:Integrations>1</nci:Integrations>
                    <nci:EtcId></nci:EtcId>
                </nci:FilterConfig>
            </nci:Filters>
        </nci:NircamImaging>
    </Template>
    <ScienceDuration>1020</ScienceDuration>
    <CoordinatedParallel>false</CoordinatedParallel>
    <SpecialRequirements>
        <Between After="{after}" Before="{before}"/>
        <OrientRange OrientMin="{eph['sunTargetPA'][i] - 3:.0f} Degrees" OrientMax="{eph['sunTargetPA'][i] + 3:.0f} Degrees" MsaSelectedAngle="false"/>
    </SpecialRequirements>
    <MosaicParameters>
        <Rows>1</Rows>
        <Columns>1</Columns>
        <RowOverlapPercent>10.0</RowOverlapPercent>
        <ColumnOverlapPercent>10.0</ColumnOverlapPercent>
        <SkewDegreesX>0.0</SkewDegreesX>
        <SkewDegreesY>0.0</SkewDegreesY>
        <MosaicTiles TileNumber="1">
            <TileState>Tile Included</TileState>
        </MosaicTiles>
        <MosaicTileOrder>DEFAULT</MosaicTileOrder>
    </MosaicParameters>
    <ObservingWindowsContainer>
        <MOSSStartDate></MOSSStartDate>
        <MOSSEndDate></MOSSEndDate>
        <MOSSShowWindows>false</MOSSShowWindows>
    </ObservingWindowsContainer>
    <Visit Number="1"/>
</Observation>
""",
                    xml_parser,
                )
            )

        if args.nirspec_ifu:
            after, before = format_between(eph[i])
            nirspec_ifu.append(
                etree.XML(
                    f"""
<Observation AutoTarget="false">
    <Number>{i + 101}</Number>
    <TargetID>{target_id}</TargetID>
    <Label>Spectrum</Label>
    <Instrument>NIRSPEC</Instrument>
    <Template xmlns:nsifus="http://www.stsci.edu/JWST/APT/Template/NirspecIFUSpectroscopy">
        <nsifus:NirspecIFUSpectroscopy>
            <nsifus:TaMethod>VERIFY_ONLY</nsifus:TaMethod>
            <nsifus:PointingVerificationImage xmlns:ns="http://www.stsci.edu/JWST/APT/Instrument/Nirspec">
                <ns:Filter>CLEAR</ns:Filter>
                <ns:ReadoutPattern>NRSIRS2RAPID</ns:ReadoutPattern>
                <ns:Groups>10</ns:Groups>
                <ns:Integrations>1</ns:Integrations>
                <ns:EtcId></ns:EtcId>
                <ns:MsaAcquisitionConfigFile>ALLCLOSED</ns:MsaAcquisitionConfigFile>
            </nsifus:PointingVerificationImage>
            <nsifus:DitherType>CYCLING</nsifus:DitherType>
            <nsifus:DitherSize>LARGE</nsifus:DitherSize>
            <nsifus:StartingPoint>1</nsifus:StartingPoint>
            <nsifus:NumberOfPoints>16</nsifus:NumberOfPoints>
            <nsifus:Exposures>
                <nsifus:Exposure>
                    <nsifus:Grating>PRISM</nsifus:Grating>
                    <nsifus:Filter>CLEAR</nsifus:Filter>
                    <nsifus:ReadoutPattern>NRSIRS2</nsifus:ReadoutPattern>
                    <nsifus:Groups>20</nsifus:Groups>
                    <nsifus:Integrations>1</nsifus:Integrations>
                    <nsifus:EtcId>59382.63</nsifus:EtcId>
                    <nsifus:Autocal>NONE</nsifus:Autocal>
                    <nsifus:Leakcal>false</nsifus:Leakcal>
                    <nsifus:Dither>true</nsifus:Dither>
                </nsifus:Exposure>
            </nsifus:Exposures>
        </nsifus:NirspecIFUSpectroscopy>
    </Template>
    <ScienceDuration>23344</ScienceDuration>
    <CoordinatedParallel>false</CoordinatedParallel>
    <SpecialRequirements>
        <Between After="{after}" Before="{before}"/>
        <OrientRange OrientMin="{eph['sunTargetPA'][i] - 3:.0f} Degrees" OrientMax="{eph['sunTargetPA'][i] + 3:.0f} Degrees" MsaSelectedAngle="false"/>
    </SpecialRequirements>
    <MosaicParameters>
        <Rows>1</Rows>
        <Columns>1</Columns>
        <RowOverlapPercent>10.0</RowOverlapPercent>
        <ColumnOverlapPercent>10.0</ColumnOverlapPercent>
        <SkewDegreesX>0.0</SkewDegreesX>
        <SkewDegreesY>0.0</SkewDegreesY>
        <MosaicTiles TileNumber="1">
            <TileState>Tile Included</TileState>
        </MosaicTiles>
        <MosaicTileOrder>DEFAULT</MosaicTileOrder>
    </MosaicParameters>
    <ObservingWindowsContainer>
        <ObservingWindowSpec xsi:type="DefaultAngularVelocityWindow" Object1="1 HALE-BOPP" Object2="" Observer="JWST" Velocity="0.03" Condition="LESS THAN" Within="Within" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"/>
        <MOSSStartDate></MOSSStartDate>
        <MOSSEndDate></MOSSEndDate>
        <MOSSShowWindows>false</MOSSShowWindows>
    </ObservingWindowsContainer>
    <Visit Number="1"/>
</Observation>
""",
                    xml_parser,
                )
            )

    if args.nircam:
        data_requests.append(nircam)

    if args.nirspec_ifu:
        data_requests.append(nirspec_ifu)

    return etree.ElementTree(root)


eph = get_ephemeris(args)

if args.xml:
    tree = eph_to_xml(eph)
    opts = dict(encoding="UTF-8", standalone=True, pretty_print=True)
    if args.o is sys.stdout:
        print(etree.tostring(tree, **opts).decode())
    else:
        tree.write(args.o, **opts)
else:
    tab = eph_to_table(eph)
    tab.write(args.o, format="ascii.commented_header", overwrite=True)
