# ==================#
# STATIONXML.PY
# Reading from Antelope database, this script creates stationXML of the stations.
# For future change keep in mind that a precise order of the attributes is required (check stationXML docs).

# ==================#
# INTRODUCTION

# ==================#
# 	TODO
#
# * Uncomment the correct argument parser (args = parser.parse_args()) and delete the others
# - Write an introduction

# ==================#
# PACKAGES

# Generic Packages
import warnings
import argparse
import xml.etree.ElementTree as ET
import configparser
from typing import List
from datetime import datetime, timezone
import re
import subprocess
from serviceModule import (
    Coefficient,
    Decimation,
    JsonEncoder,
    Level,
    PoleZero,
    PolesZeros,
    StationInfo,
)
from serviceModule import iso8601DateTime, seedIdentifier, latitudeFloat, longitudeFloat
from serviceModule import filterTable

warnings.filterwarnings("ignore")

# Antelope packages
from antelope.datascope import *

# ==================#
# DEFINITIONS
config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(sys.argv[0]), "config.ini"))
DB_PATH = config["DEFAULT"]["DB_PATH"]
DB_NAME = config["DEFAULT"]["DB_NAME"]
PROGRAM_VERSION = config["PROGRAM_INFO"]["PROGRAM_VERSION"]
VALIDATOR = config["DEFAULT"]["VALIDATOR_PATH"]

NM2M = True  # Flag to perform the conversion from nanometer to meter
SKIPCHAN = ["CPU", "L.."]  # Channel name pattern to skip


def writeRawStation(
    stationInfos: List[StationInfo], limit: int, skip: int, output: str = None
) -> None:
    index = 0
    printed = 0

    # with open("stationxml.json", "w") as f:
    #  	f.write(JsonEncoder().encode(stationInfos))

    # stationJson = JsonEncoder().encode(stationInfos)
    # stationInfoDecoded = json.loads(stationJson, object_hook=lambda d: SimpleNamespace(**d))
    # c = 5

    if output:
        with open(output, "w") as f:
            for stationInfo in stationInfos:
                index = index + 1
                if index <= skip:
                    continue

                if printed >= limit:
                    return

                f.write(JsonEncoder().encode(stationInfo))
                printed = printed + 1
    else:
        for stationInfo in stationInfos:
            index = index + 1
            if index <= skip:
                continue

            if printed >= limit:
                return

            print(JsonEncoder().encode(stationInfo))
            printed = printed + 1


def _pretty_print(current, parent=None, index=-1, depth=0):
    """Function to format the XML output."""
    for i, node in enumerate(current):
        _pretty_print(node, current, i, depth + 1)
    if parent is not None:
        if index == 0:
            parent.text = "\n" + ("\t" * depth)
        else:
            parent[index - 1].tail = "\n" + ("\t" * depth)
        if index == len(parent) - 1:
            current.tail = "\n" + ("\t" * (depth - 1))


# ==================#
# ARGUMENTS

parser = argparse.ArgumentParser(
    description="Writes in the terminal (or into a file) an XML, conformed to version 1.2 of StationXML FDSN specification",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

# Temporal Constraints
startTime = None
endTime = None
startBefore = None
startAfter = None
endBefore = None
endAfter = None
parser.add_argument(
    "--starttime",
    type=iso8601DateTime,
    help="Limit metadata to epochs starting on or after this time",
)
parser.add_argument(
    "--endtime", type=iso8601DateTime, help="Limit metadata to epochs ending on or before this time"
)
parser.add_argument(
    "--startbefore",
    type=iso8601DateTime,
    help="Limit metadata to stations starting on or before this time",
)
parser.add_argument(
    "--startafter",
    type=iso8601DateTime,
    help="Limit metadata to stations starting on or after this time",
)
parser.add_argument(
    "--endbefore",
    type=iso8601DateTime,
    help="Limit metadata to stations ending on or before this time",
)
parser.add_argument(
    "--endafter",
    type=iso8601DateTime,
    help="Limit metadata to stations ending on or after this time",
)

# SEED Identifiers
networks = []
stations = []
locations = []
channels = []
parser.add_argument(
    "--network",
    type=seedIdentifier,
    help="SEED Network code. Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)",
)
parser.add_argument(
    "--station",
    type=seedIdentifier,
    help="SEED Station code. Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)",
)
parser.add_argument(
    "--location",
    type=seedIdentifier,
    help="SEED Location code (use '--' for blank location identifiers). Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)",
)
parser.add_argument(
    "--channel",
    type=seedIdentifier,
    help="SEED Channel code. Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)",
)

# Rectangular Spatial Selection
minLatitude = -90.0
maxLatitude = 90.0
minLongitude = -180.0
maxLongitude = 180.0
parser.add_argument("--minlatitude", type=latitudeFloat, help="Southern boundary", default=-90.0)
parser.add_argument("--maxlatitude", type=latitudeFloat, help="Northern boundary", default=90.0)
parser.add_argument("--minlongitude", type=longitudeFloat, help="Western boundary", default=-180.0)
parser.add_argument("--maxlongitude", type=longitudeFloat, help="Eastern boundary", default=180.0)

# Circular Spatial Selection
latitude = 0.0
longitude = 0.0
maxRadius = 0.0
minRadius = 0.0
parser.add_argument("--latitude", type=latitudeFloat, help="Central latitude point", default=0.0)
parser.add_argument("--longitude", type=longitudeFloat, help="Central longitude point", default=0.0)
parser.add_argument(
    "--maxradius",
    type=float,
    help="Maximum radial arc-distance from the central point	",
    default=0.0,
)
parser.add_argument(
    "--minradius",
    type=float,
    help="Minimum radial arc-distance from the central point",
    default=0.0,
)

# Request Options
level = "station"
includedRestricted = True
parser.add_argument(
    "--level",
    type=str,
    help="Level of detail ('network', 'station', 'channel', or 'response')",
    default="station",
    choices=["network", "station", "channel", "response"],
)
parser.add_argument(
    "--includerestricted",
    type=bool,
    help="Include information on restricted stations",
    default=True,
)

# Version
parser.add_argument(
    "-V",
    "--version",
    help="Show the version of the program and exit",
    action="version",
    version=PROGRAM_VERSION,
)
# parser.add_argument('-V', '--version', help='Show the version of the program and exit', action='version', version='%(prog)s {0}'.format(PROGRAM_VERSION))

# Other
raw = False
limit = 10000
skip = 0
output = ""
parser.add_argument(
    "-R",
    "--raw",
    help="Prints the raw data, after being filtered, in a JSON format",
    action="store_true",
)
parser.add_argument(
    "-L",
    "--limit",
    type=int,
    help="Limits the number of results. Used together with --raw. Default is 10000",
    default=10000,
)
parser.add_argument(
    "-S",
    "--skip",
    type=int,
    help="Skips number of results. Used together with --raw. Default is 0",
    default=0,
)
parser.add_argument("-O", "--output", type=str, help="Saves the result in a file with this name")

# Parse the arguments
args = parser.parse_args()
# args = parser.parse_args(['--network', 'RF,I?', '--starttime', '2003-03-12'])
# args = parser.parse_args(['--station', 'CESC', '--level', 'station', '--raw'])

if args.starttime:
    startTime = args.starttime
if args.endtime:
    endTime = args.endtime
if args.startbefore:
    startBefore = args.startbefore
if args.startafter:
    startAfter = args.startafter
if args.endbefore:
    endBefore = args.endbefore
if args.endafter:
    endAfter = args.endafter

if args.network:
    networks = args.network
if args.station:
    stations = args.station
if args.location:
    locations = []
    # The value '--' is used when the location is blank. So we
    # must remove all entries with that value and leave them blank too
    # This doesn't happen for RegExes
    for location in args.location:
        if location[0] == False:
            # '--' can't be normally passed directly, as it is otherwise interpreted as a broken argument
            # So we must check any combination of quotes
            if location[1] == "--" or location[1] == "'--'" or location[1] == '"--"':
                locations.append((False, ""))
            else:
                # Not a blank location
                locations.append(location)
        else:
            # Is a RegEx
            locations.append(location)
if args.channel:
    channels = args.channel

if args.minlatitude:
    minLatitude = args.minlatitude
if args.maxlatitude:
    maxLatitude = args.maxlatitude
if args.minlongitude:
    minLongitude = args.minlongitude
if args.maxlongitude:
    maxLongitude = args.maxlongitude

if args.latitude:
    latitude = args.latitude
if args.longitude:
    longitude = args.longitude
if args.maxradius:
    maxRadius = args.maxradius
if args.minradius:
    minRadius = args.minradius

if args.level:
    level = args.level
if args.includerestricted:
    # Unused in the code
    includedRestricted = args.includerestricted

if args.raw:
    raw = args.raw
if args.limit:
    limit = args.limit
if args.skip:
    skip = args.skip
if args.output:
    output = args.output


# ==================#
# IMPLEMENTATION

REGCHAN = [re.compile(uc) for uc in SKIPCHAN]

# DB list in server (eg. 140.105.54.8) on path (eg. /home/rt/rtsystem/db/)
# Open database in read mode
db = dbopen(DB_PATH + DB_NAME, perm="r")

# print(DB_PATH + DB_NAME)

# Prepare the table structure
# Filter the results as-you-go, to reduce the amount of data
# between each step
networkTable = db.lookup(table="network")
# print(len(networks))

if len(networks) > 0:
    networkTable = filterTable(networkTable, "network.net", networks)
    # Per specification, networks should be filtrable by start/end date.
    # None is present in the database.

# print(networkTable.find_join_tables())
joinAffiliationTable = networkTable.process(["dbjoin affiliation"]).subset(
    "affiliation.net == network.net"
)
if len(stations) > 0:
    joinAffiliationTable = filterTable(joinAffiliationTable, "affiliation.sta", stations)

joinSiteTable = joinAffiliationTable.process(["dbjoin site"]).subset("site.sta == affiliation.sta")
if Level[level] == Level.station:
    # As per specification, filter by date only at the specified level
    if startTime:
        joinSiteTable = joinSiteTable.subset(
            "site.ondate >= {0}".format(startTime.strftime("%Y%j"))
        )
    if endTime:
        joinSiteTable = joinSiteTable.subset("site.offdate <= {0}".format(endTime.strftime("%Y%j")))
    if startBefore:
        joinSiteTable = joinSiteTable.subset(
            "site.ondate < {0}".format(startBefore.strftime("%Y%j"))
        )
    if startAfter:
        joinSiteTable = joinSiteTable.subset(
            "site.ondate > {0}".format(startAfter.strftime("%Y%j"))
        )
    if endBefore:
        # In the table, stations without and ending date have offdate set to -1. These must be excluded
        joinSiteTable = joinSiteTable.subset(
            "site.offdate < {0} && site.offdate > 0".format(endBefore.strftime("%Y%j"))
        )
    if endAfter:
        # In the table, stations without and ending date have offdate set to -1. These must be included
        joinSiteTable = joinSiteTable.subset(
            "site.offdate > {0} || site.offdate < 0".format(endAfter.strftime("%Y%j"))
        )

if Level[level].value > Level.network.value:
    # Note that the specification defines that this filter should be applied according to the
    # level, but the only reference to spatial coordinates in the databse is found in the station.

    # Default values for Rectangular Spatial Selection cover the entire globe,
    # so perform a filter only if they are different from their defaults.
    if minLatitude != -90.0:
        joinSiteTable = joinSiteTable.subset("site.lat >= {0}".format(minLatitude))
    if maxLatitude != 90.0:
        joinSiteTable = joinSiteTable.subset("site.lat <= {0}".format(maxLatitude))
    if minLongitude != -180.0:
        joinSiteTable = joinSiteTable.subset("site.lon >= {0}".format(minLongitude))
    if maxLongitude != 180.0:
        joinSiteTable = joinSiteTable.subset("site.lon <= {0}".format(maxLongitude))

    # Default values for Circula Spatial Selection filter out everything,
    # so perform a filter only if they are different from their defaults
    if latitude > 0.0 or longitude > 0.0 or maxRadius > 0.0 or minRadius > 0.0:
        # The distance is calculated using the Pythagorean theorem.
        # The center is located in (latitude, longitude).
        # The position of each station relative to the center is then (site.lat - latitude, site.lon - longitude).
        # The Eucliean distance between the center and the position of the station must then fall between
        # minRadius and maxRadius.
        joinSiteTable = joinSiteTable.subset(
            "(((site.lat - {0}) * (site.lat - {0})) + ((site.lon - {1}) * (site.lon - {1}))) <= {2} * {2} && (((site.lat - {0}) * (site.lat - {0})) + ((site.lon - {1}) * (site.lon - {1}))) >= {3} * {3}".format(
                latitude, longitude, maxRadius, minRadius
            )
        )

joinSChanLocTable = joinSiteTable.process(["dbjoin schanloc"]).subset("site.sta == schanloc.sta")
if len(channels) > 0:
    joinSChanLocTable = filterTable(joinSChanLocTable, "schanloc.fchan", channels)
if len(locations) > 0:
    joinSChanLocTable = filterTable(joinSChanLocTable, "schanloc.loc", locations)

joinSiteChanTable = joinSChanLocTable.process(["dbjoin sitechan"]).subset(
    "sitechan.chan == schanloc.chan"
)
if Level[level].value > Level.station.value:
    # As per specification, filter by date only at the specified level
    if startTime:
        joinSiteChanTable = joinSiteChanTable.subset(
            "sitechan.ondate >= {0}".format(startTime.strftime("%Y%j"))
        )
    if endTime:
        joinSiteChanTable = joinSiteChanTable.subset(
            "sitechan.offdate <= {0}".format(endTime.strftime("%Y%j"))
        )
    if startBefore:
        joinSiteChanTable = joinSiteChanTable.subset(
            "sitechan.ondate < {0}".format(startBefore.strftime("%Y%j"))
        )
    if startAfter:
        joinSiteChanTable = joinSiteChanTable.subset(
            "sitechan.ondate > {0}".format(startAfter.strftime("%Y%j"))
        )
    if endBefore:
        # In the table, channels without and ending date have offdate set to -1. These must be excluded
        joinSiteChanTable = joinSiteChanTable.subset(
            "sitechan.offdate < {0} && sitechan.offdate > 0".format(endBefore.strftime("%Y%j"))
        )
    if endAfter:
        # In the table, channels without and ending date have offdate set to -1. These must be included
        joinSiteChanTable = joinSiteChanTable.subset(
            "sitechan.offdate > {0} || sitechan.offdate < 0".format(endAfter.strftime("%Y%j"))
        )

joinSensor = joinSiteChanTable.process(["dbjoin sensor"]).subset("sensor.chan == schanloc.chan")
table = joinSensor.process(["dbjoin calibration"]).subset("calibration.chan == schanloc.chan")
if Level[level] == Level.response or raw:
    # Join Stage table only if the level is 'response' or if raw data is requested. This is done only because
    # this table isn't filtered by any SEED identifier but appears only if the level is high enough.
    # Done only for performance reasons. Feel free to remove the 'if' statement if you think this
    # won't impact too much the processing time
    table = table.process(["dbjoin stage"]).subset("stage.chan == schanloc.chan")
    # Sort the table
    table = table.sort(["network.net", "stage.sta", "stage.chan", "stage.time", "stage.stageid"])
else:
    # Sort the table with the available fields
    table = table.sort(["network.net", "site.sta", "sitechan.chan"])

previousNetwork = ""  # {Network code}
previousStation = ""  # {Station Code}_{Station Start Date}
previousChannel = ""  # {Station Code}_{Channel+Location Code}_{Channel+Location Start Date}
previousResponse = ""  # {Station Code}_{Channel+Location Code}_{Calibration Start Epoch}
previousStage = ""  # {Station Code}_{Channel+Location Code}_{Calibration Start Epoch}_{Stage Id}
previousXmlNetwork = None
previousXmlStation = None
previousXmlChannel = None
previousXmlResponse = None
previousXmlStage = None

# Check the number of elements. If there's none,
# return nothing.
nrec = table.query(dbRECORD_COUNT)
# print(nrec)
if nrec == 0:
    if len(output) > 0:
        # Create an empty file
        open(output, "a").close()
    else:
        print()

    exit()

# Create XML header
xmlRoot = ET.Element("FDSNStationXML")
xmlRoot.set("xmlns", "http://www.fdsn.org/xml/station/1")
xmlRoot.set("schemaVersion", "1.2")
ET.SubElement(xmlRoot, "Source").text = "RF"
ET.SubElement(xmlRoot, "Sender").text = "University of Trieste - SeisRaM"  # SFF
ET.SubElement(xmlRoot, "Module").text = "stationxml.py"
# ET.SubElement(xmlRoot, 'ModuleURI').text = 'PLACEHOLDER'
ET.SubElement(xmlRoot, "Created").text = datetime.utcnow().strftime("%Y-%m-%dT%H:%m:%S.%fZ")

if raw:
    stationInfos = []

    for i in range(nrec):
        table.record = i
        (
            networkCode,
            networkName,
            stationCode,
            stationName,
            stationLatitude,
            stationLongitude,
            stationElevation,
            stationStartDate,
            stationEndDate,
            channelCode,
            channelForeignCode,
            locationCode,
            channelStartDate,
            channelEndDate,
            channelDepth,
            channelDescription,
            calibrationInstrumentName,
            calibrationSegType,
            calibrationSampleRate,
            calibrationChannel,
            calibrationFrequency,
            calibrationUnits,
            calibrationTime,
            calibrationEndTime,
            calibrationCalibration,
            calibrationCalper,
            stageChannel,
            stageTime,
            stageEndTime,
            stageStageId,
            stageGnom,
            stageInputUnits,
            stageOutputUnits,
            stageSampleRate,
            stageGainCalibration,
            stageDecifac,
            stageIzero,
        ) = table.getv(
            "network.net",
            "network.netname",
            "site.sta",
            "site.staname",
            "site.lat",
            "site.lon",
            "site.elev",
            "site.ondate",
            "site.offdate",
            "schanloc.chan",
            "schanloc.fchan",
            "schanloc.loc",
            "sitechan.ondate",
            "sitechan.offdate",
            "sitechan.edepth",
            "sitechan.descrip",
            "calibration.insname",
            "calibration.segtype",
            "calibration.samprate",
            "calibration.chan",
            "calibration.fc",
            "calibration.units",
            "calibration.time",
            "calibration.endtime",
            "calibration.calib",
            "calibration.calper",
            "stage.chan",
            "stage.time",
            "stage.endtime",
            "stage.stageid",
            "stage.gnom",
            "stage.iunits",
            "stage.ounits",
            "stage.samprate",
            "stage.gcalib",
            "stage.decifac",
            "stage.izero",
        )
        poles = []
        zeros = []
        coefficients = []
        polesZeros = None
        decimation = None

        filename = table.filename()
        if filename[0] == 1:
            with open(filename[1], "r") as file:
                # Remove header comments
                metadata = ""
                while True:
                    line = file.readline()
                    if not line.startswith("#"):
                        metadata = line.split()
                        break

                if metadata[3] == "paz":
                    # Poles and Zeros

                    normalizationFactor = (
                        file.readline().replace("\n", "").replace("\r", "").split()[0]
                    )
                    mode = ""
                    number = 0
                    while True:
                        line = file.readline()
                        if len(line) == 0:
                            break

                        splitLine = line.split()
                        if len(splitLine) == 2:
                            mode = splitLine[1]
                            number = 0
                        elif len(splitLine) == 4:
                            if mode == "Poles":
                                poles.append(
                                    PoleZero(
                                        number,
                                        float(splitLine[0]),
                                        float(splitLine[2]),
                                        float(splitLine[1]),
                                        float(splitLine[3]),
                                    )
                                )
                                number = number + 1
                            elif mode == "Zeros":
                                zeros.append(
                                    PoleZero(
                                        number,
                                        float(splitLine[0]),
                                        float(splitLine[2]),
                                        float(splitLine[1]),
                                        float(splitLine[3]),
                                    )
                                )
                                number = number + 1

                    polesZeros = PolesZeros(normalizationFactor, poles, zeros)
                elif metadata[3] == "fir":
                    # FIR
                    decimationString = file.readline().split()
                    decimation = Decimation(decimationString[0], decimationString[1])

                    while True:
                        line = file.readline()
                        if len(line) == 0:
                            break

                        splitLine = line.split()
                        if len(splitLine) == 1:
                            continue
                        elif len(splitLine) == 2:
                            coefficients.append(
                                Coefficient(float(splitLine[0]), float(splitLine[1]))
                            )

        else:
            decimation = Decimation(str(stageSampleRate), stageDecifac)

        stationInfos.append(
            StationInfo(
                networkCode,
                networkName,
                stationCode,
                stationName,
                stationLatitude,
                stationLongitude,
                stationElevation,
                stationStartDate,
                stationEndDate,
                channelCode,
                channelForeignCode,
                locationCode,
                channelStartDate,
                channelEndDate,
                channelDepth,
                channelDescription,
                calibrationInstrumentName,
                calibrationSegType,
                calibrationSampleRate,
                calibrationChannel,
                calibrationFrequency,
                calibrationUnits,
                calibrationTime,
                calibrationEndTime,
                calibrationCalibration,
                calibrationCalper,
                stageChannel,
                stageTime,
                stageEndTime,
                stageStageId,
                stageGnom,
                stageInputUnits,
                stageOutputUnits,
                stageSampleRate,
                stageGainCalibration,
                stageDecifac,
                stageIzero,
                polesZeros,
                decimation,
                coefficients,
            )
        )
    writeRawStation(stationInfos, limit, skip, output)

else:
    for i in range(nrec):
        table.record = i

        # Create Network level XML Tags
        networkInfo = table.getv("network.net", "network.netname")
        if previousNetwork != networkInfo[0]:
            xmlNetwork = ET.SubElement(xmlRoot, "Network")
            xmlNetwork.set("code", networkInfo[0])
            ET.SubElement(xmlNetwork, "Description").text = networkInfo[1]

            previousNetwork = networkInfo[0]
            previousXmlNetwork = xmlNetwork

        if Level[level].value < Level.station.value:
            continue

        # Create Station level XML Tags
        stationInfo = table.getv(
            "site.sta",
            "site.staname",
            "site.lat",
            "site.lon",
            "site.elev",
            "site.ondate",
            "site.offdate",
        )
        if previousStation != f"{stationInfo[0]}_{stationInfo[5]}":
            xmlStation = ET.SubElement(previousXmlNetwork, "Station")
            xmlStation.set("code", stationInfo[0])
            # xmlStation.set("startDate", str(datetime.strptime(str(stationInfo[5]), "%Y%j"))) #OG
            xmlStation.set(
                "startDate",
                datetime.strptime(str(stationInfo[5]), "%Y%j").strftime("%Y-%m-%dT%H:%m:%S.%fZ"),
            )  # SFF
            if stationInfo[6] > 0:
                # xmlStation.set("endDate", str(datetime.strptime(str(stationInfo[6]), "%Y%j"))) # OG
                xmlStation.set(
                    "endtDate",
                    datetime.strptime(str(stationInfo[6]), "%Y%j").strftime(
                        "%Y-%m-%dT%H:%m:%S.%fZ"
                    ),
                )  # SFF
            ET.SubElement(xmlStation, "Latitude").text = str(stationInfo[2])
            ET.SubElement(xmlStation, "Longitude").text = str(stationInfo[3])
            ET.SubElement(xmlStation, "Elevation").text = str(
                stationInfo[4] * 1000
            )  # Elevation in Antelope is in Km
            # ET.SubElement(xmlStation, 'Vault').text = 'PLACEHOLDER'
            # ET.SubElement(xmlStation, "CreationDate").text = str(
            # datetime.strptime(str(stationInfo[5]), "%Y%j")
            # )
            # ET.SubElement(xmlStation, "CreationDate").text = datetime.strptime(
            #     str(stationInfo[5]), "%Y%j"
            # ).strftime(
            #     "%Y-%m-%dT%H:%m:%S.%fZ"
            # )  # SFF

            xmlStationSite = ET.SubElement(xmlStation, "Site")
            ET.SubElement(xmlStationSite, "Name").text = stationInfo[1]

            previousStation = f"{stationInfo[0]}_{stationInfo[5]}"
            previousXmlStation = xmlStation

        if Level[level].value < Level.channel.value:
            continue

        # Create Channel level XML Tags
        channelInfo = table.getv(
            "schanloc.chan",
            "schanloc.fchan",
            "schanloc.loc",
            "site.lat",
            "site.lon",
            "site.elev",
            "sitechan.ondate",
            "sitechan.offdate",
            "sitechan.edepth",
            "sitechan.vang",
            "sitechan.hang",
            "sitechan.descrip",
            "calibration.insname",
            "calibration.segtype",
            "calibration.samprate",
        )
        if any(bool(r.match(channelInfo[0])) for r in REGCHAN):
            continue
        if previousChannel != f"{stationInfo[0]}_{channelInfo[0]}_{channelInfo[6]}":
            xmlChannel = ET.SubElement(previousXmlStation, "Channel")
            xmlChannel.set("code", channelInfo[1])
            xmlChannel.set("locationCode", channelInfo[2])
            # xmlChannel.set("startDate", str(datetime.strptime(str(channelInfo[6]), "%Y%j"))) #OG
            xmlChannel.set(
                "startDate",
                datetime.strptime(str(channelInfo[6]), "%Y%j").strftime("%Y-%m-%dT%H:%m:%S.%fZ"),
            )  # SFF
            if channelInfo[7] > 0:
                # xmlChannel.set("endDate", str(datetime.strptime(str(channelInfo[7]), "%Y%j"))) #OG
                xmlChannel.set(
                    "endDate",
                    datetime.strptime(str(channelInfo[7]), "%Y%j").strftime(
                        "%Y-%m-%dT%H:%m:%S.%fZ"
                    ),
                )  # SFF
            xmlChannelComment = ET.SubElement(xmlChannel, "Comment")
            ET.SubElement(xmlChannelComment, "Value").text = channelInfo[11]

            ET.SubElement(xmlChannel, "Latitude").text = str(channelInfo[3])
            ET.SubElement(xmlChannel, "Longitude").text = str(channelInfo[4])
            ET.SubElement(xmlChannel, "Elevation").text = str(
                channelInfo[5] * 1000
            )  # Elevation in Antelope is in Km
            ET.SubElement(xmlChannel, "Depth").text = str(channelInfo[8])
            if channelInfo[10] != -999.9:
                ET.SubElement(xmlChannel, "Azimuth").text = str(
                    channelInfo[10]
                )  # TODO Check if it is correct
            # In Antelope: vhang is 90, 0, 180 for horizontal, upward vertical, and downward vertical component, respectively
            # In stationXML: dip is 0, -90, 90 for horizontal, upward vertical, and downward vertical component, respectively
            if channelInfo[9] != -999.9:
                ET.SubElement(xmlChannel, "Dip").text = str(
                    channelInfo[9] - 90
                )  # TODO Check if it is correct
            ET.SubElement(xmlChannel, "SampleRate").text = str(channelInfo[14])

            xmlChannelSensor = ET.SubElement(xmlChannel, "Sensor")
            ET.SubElement(xmlChannelSensor, "Type").text = channelInfo[13]
            ET.SubElement(xmlChannelSensor, "Description").text = channelInfo[12]

            previousChannel = f"{stationInfo[0]}_{channelInfo[0]}_{channelInfo[6]}"
            previousXmlChannel = xmlChannel

        if Level[level].value < Level.response.value:
            continue

        # Create Response level XML Tags
        responseInfo = table.getv(
            "calibration.chan",
            "calibration.fc",
            "calibration.units",
            "calibration.time",
            "calibration.endtime",
            "calibration.calib",
        )

        # Ad hoc conversion from nm (or nm/s or nm/s**2) to m (or m/s or m/s**2)
        if NM2M and "nm" in responseInfo[2]:
            tmpresponseInfo = list(responseInfo)
            tmpresponseInfo[2] = tmpresponseInfo[2].replace("nm", "m")
            tmpresponseInfo[5] *= 1e-9
            responseInfo = tuple(tmpresponseInfo)

        if previousResponse != f"{stationInfo[0]}_{responseInfo[0]}_{responseInfo[3]}":
            xmlResponse = ET.SubElement(previousXmlChannel, "Response")

            xmlResponseInstrumentSensitivity = ET.SubElement(xmlResponse, "InstrumentSensitivity")
            ET.SubElement(xmlResponseInstrumentSensitivity, "Value").text = str(
                1 / float(responseInfo[5])
            )
            ET.SubElement(xmlResponseInstrumentSensitivity, "Frequency").text = str(responseInfo[1])

            xmlResponseInputUnits = ET.SubElement(xmlResponseInstrumentSensitivity, "InputUnits")
            ET.SubElement(xmlResponseInputUnits, "Name").text = responseInfo[2]

            xmlResponseOutputUnits = ET.SubElement(xmlResponseInstrumentSensitivity, "OutputUnits")
            ET.SubElement(xmlResponseOutputUnits, "Name").text = (
                "count"  # TODO Not sure. There is no output units.
            )

            previousResponse = f"{stationInfo[0]}_{responseInfo[0]}_{responseInfo[3]}"
            previousXmlResponse = xmlResponse
        # Create Stage level XML Tags. It is part of the Response
        stageInfo = table.getv(
            "stage.chan",
            "stage.time",
            "stage.endtime",
            "stage.stageid",
            "stage.gnom",
            "calibration.fc",
            "stage.iunits",
            "stage.ounits",
            "stage.samprate",
            "stage.gcalib",
            "calibration.fc",
            "calibration.calper",
            "stage.decifac",
            "stage.izero",
        )

        # Ad hoc conversion from nm (or nm/s or nm/s**2) to m (or m/s or m/s**2)
        if NM2M and "nm" in stageInfo[6]:
            tmpstageInfo = list(stageInfo)
            tmpstageInfo[4] *= 1e9
            tmpstageInfo[6] = tmpstageInfo[6].replace("nm", "m")
            stageInfo = tuple(tmpstageInfo)

        if previousStage != f"{stationInfo[0]}_{stageInfo[0]}_{stageInfo[1]}_{stageInfo[3]}":
            xmlStage = ET.SubElement(previousXmlResponse, "Stage")
            xmlStage.set("number", str(stageInfo[3]))

            # Get the linked file and parse it
            filename = table.filename()
            if filename[0] == 1:
                with open(filename[1], "r") as file:
                    # Remove header comments
                    metadata = ""
                    while True:
                        line = file.readline()
                        if not line.startswith("#"):
                            metadata = line.split()
                            break

                    if metadata[3] == "paz":
                        # Poles and Zeros
                        xmlPolesZeros = ET.SubElement(xmlStage, "PolesZeros")

                        xmlPolesZerosInputUnits = ET.SubElement(xmlPolesZeros, "InputUnits")
                        inputName = stageInfo[6]
                        if len(inputName) == 0 or inputName == "-":
                            inputName = "count"  # Hardcoded as 'count' if none is present
                        ET.SubElement(xmlPolesZerosInputUnits, "Name").text = inputName

                        xmlPolesZerosOutputUnits = ET.SubElement(xmlPolesZeros, "OutputUnits")
                        outputName = stageInfo[7]
                        if len(outputName) == 0 or outputName == "-" or outputName == "counts":
                            outputName = "count"  # Hardcoded as 'count' if none is present
                        ET.SubElement(xmlPolesZerosOutputUnits, "Name").text = outputName

                        ET.SubElement(xmlPolesZeros, "PzTransferFunctionType").text = (
                            "LAPLACE (RADIANS/SECOND)"  # Hardcoded
                        )
                        normalizationFactor = (
                            file.readline().replace("\n", "").replace("\r", "").split()[0]
                        )
                        ET.SubElement(xmlPolesZeros, "NormalizationFactor").text = (
                            normalizationFactor
                        )
                        ET.SubElement(xmlPolesZeros, "NormalizationFrequency").text = str(
                            stageInfo[11]
                        )

                        mode = ""
                        number = 0
                        while True:
                            line = file.readline()
                            if len(line) == 0:
                                break

                            splitLine = line.split()
                            if len(splitLine) == 2:
                                mode = splitLine[1]
                                number = 0
                            elif len(splitLine) == 4:
                                if mode == "Zeros":
                                    xmlZero = ET.SubElement(xmlPolesZeros, "Zero")
                                    xmlZero.set("number", str(number))

                                    xmlZeroReal = ET.SubElement(xmlZero, "Real")
                                    xmlZeroReal.text = str(float(splitLine[0]))

                                    error = float(splitLine[2])
                                    xmlZeroReal.set("minusError", str(min([0.0, error])))
                                    xmlZeroReal.set("plusError", str(max([0.0, error])))
                                    # if error < 0:
                                    #     xmlZeroReal.set("plusError", str(0.0))
                                    #     xmlZeroReal.set("minusError", str(error))
                                    # else:
                                    #     xmlZeroReal.set("minusError", str(0.0))
                                    #     xmlZeroReal.set("plusError", str(error))

                                    xmlZeroImaginary = ET.SubElement(xmlZero, "Imaginary")
                                    xmlZeroImaginary.text = str(float(splitLine[1]))

                                    error = float(splitLine[3])
                                    xmlZeroImaginary.set("minusError", str(min([0.0, error])))
                                    xmlZeroImaginary.set("plusError", str(max([0.0, error])))
                                    # if error < 0:
                                    #     xmlZeroImaginary.set("minusError", str(error))
                                    #     xmlZeroImaginary.set("plusError", str(0.0))
                                    # else:
                                    #     xmlZeroImaginary.set("minusError", str(0.0))
                                    #     xmlZeroImaginary.set("plusError", str(error))
                                    number = number + 1

                                elif mode == "Poles":
                                    xmlPole = ET.SubElement(xmlPolesZeros, "Pole")
                                    xmlPole.set("number", str(number))

                                    xmlPoleReal = ET.SubElement(xmlPole, "Real")
                                    xmlPoleReal.text = str(float(splitLine[0]))

                                    error = float(splitLine[2])
                                    xmlPoleReal.set("minusError", str(min([0.0, error])))
                                    xmlPoleReal.set("plusError", str(max([0.0, error])))
                                    # if error < 0:
                                    #     xmlPoleReal.set("plusError", str(0.0))
                                    #     xmlPoleReal.set("minusError", str(error))
                                    # else:
                                    #     xmlPoleReal.set("minusError", str(0.0))
                                    #     xmlPoleReal.set("plusError", str(error))

                                    xmlPoleImaginary = ET.SubElement(xmlPole, "Imaginary")
                                    xmlPoleImaginary.text = str(float(splitLine[1]))

                                    error = float(splitLine[3])
                                    xmlPoleImaginary.set("minusError", str(min([0.0, error])))
                                    xmlPoleImaginary.set("plusError", str(max([0.0, error])))
                                    # if error < 0:
                                    #     xmlPoleImaginary.set("minusError", str(error))
                                    #     xmlPoleImaginary.set("plusError", str(0.0))
                                    # else:
                                    #     xmlPoleImaginary.set("minusError", str(0.0))
                                    #     xmlPoleImaginary.set("plusError", str(error))
                                    number = number + 1
                        xmlPolesZeros[:] = (
                            [sub for sub in xmlPolesZeros if sub.tag not in ["Pole", "Zero"]]
                            + [sub for sub in xmlPolesZeros if sub.tag == "Zero"]
                            + [sub for sub in xmlPolesZeros if sub.tag == "Pole"]
                        )

                    elif metadata[3] == "fir":
                        # FIR
                        decimation = file.readline().split()

                        xmlCoefficients = ET.SubElement(xmlStage, "Coefficients")
                        xmlCoefficientsInputUnits = ET.SubElement(xmlCoefficients, "InputUnits")
                        inputName = stageInfo[6]
                        if len(inputName) == 0 or inputName == "-":
                            inputName = "count"  # Hardcoded as 'count' if none is present
                        ET.SubElement(xmlCoefficientsInputUnits, "Name").text = inputName
                        xmlCoefficientsOutputUnits = ET.SubElement(xmlCoefficients, "OutputUnits")
                        outputName = stageInfo[7]
                        if len(outputName) == 0 or outputName == "-" or outputName == "counts":
                            outputName = "count"  # Hardcoded as 'count' if none is present
                        ET.SubElement(xmlCoefficientsOutputUnits, "Name").text = outputName
                        ET.SubElement(xmlCoefficients, "CfTransferFunctionType").text = (
                            "DIGITAL"  # Hardcoded
                        )
                        while True:
                            line = file.readline()
                            if len(line) == 0:
                                break
                            splitLine = line.split()
                            if len(splitLine) == 1:
                                continue
                            elif len(splitLine) == 2:
                                xmlCoefficientsNumerator = ET.SubElement(
                                    xmlCoefficients, "Numerator"
                                )
                                xmlCoefficientsNumerator.text = str(float(splitLine[0]))

                                error = float(splitLine[1])
                                xmlCoefficientsNumerator.set("minusError", str(min([0.0, error])))
                                xmlCoefficientsNumerator.set("plusError", str(max([0.0, error])))
                                # if error < 0:
                                #     xmlCoefficientsNumerator.set("minusError", str(error))
                                #     xmlCoefficientsNumerator.set("plusError", str(0.0))
                                # else:
                                #     xmlCoefficientsNumerator.set("minusError", str(0.0))
                                #     xmlCoefficientsNumerator.set("plusError", str(error))

                        xmlDecimation = ET.SubElement(xmlStage, "Decimation")
                        ET.SubElement(xmlDecimation, "InputSampleRate").text = decimation[0]
                        ET.SubElement(xmlDecimation, "Factor").text = decimation[1]
                        ET.SubElement(xmlDecimation, "Offset").text = str(stageInfo[13])
                        # TODO Not sure
                        # delay = float(decimation[2]) / float(decimation[0])
                        # ET.SubElement(xmlDecimation, 'Delay').text = str(delay)		# TODO Not sure
                        # ET.SubElement(xmlDecimation, 'Correction').text = str(delay)	# TODO Not sure
                        ET.SubElement(xmlDecimation, "Delay").text = str(stageInfo[9])
                        ET.SubElement(xmlDecimation, "Correction").text = str(stageInfo[9])
            else:

                xmlCoefficients = ET.SubElement(xmlStage, "Coefficients")
                xmlCoefficientsInputUnits = ET.SubElement(xmlCoefficients, "InputUnits")
                inputName = stageInfo[6]
                if len(inputName) == 0 or inputName == "-":
                    inputName = "count"  # Hardcoded as 'count' if none is present
                ET.SubElement(xmlCoefficientsInputUnits, "Name").text = inputName
                xmlCoefficientsOutputUnits = ET.SubElement(xmlCoefficients, "OutputUnits")
                outputName = stageInfo[7]
                if len(outputName) == 0 or outputName == "-" or outputName == "counts":
                    outputName = "count"  # Hardcoded as 'count' if none is present
                ET.SubElement(xmlCoefficientsOutputUnits, "Name").text = outputName
                ET.SubElement(xmlCoefficients, "CfTransferFunctionType").text = (
                    "DIGITAL"  # Hardcoded
                )

                xmlDecimation = ET.SubElement(xmlStage, "Decimation")
                if stageInfo[8] == -1:
                    ET.SubElement(xmlDecimation, "InputSampleRate").text = str(1)
                else:
                    ET.SubElement(xmlDecimation, "InputSampleRate").text = str(stageInfo[8])
                if stageInfo[12] == -1:
                    ET.SubElement(xmlDecimation, "Factor").text = str(1)
                else:
                    ET.SubElement(xmlDecimation, "Factor").text = str(int(stageInfo[12]))
                ET.SubElement(xmlDecimation, "Offset").text = str(stageInfo[13])
                ET.SubElement(xmlDecimation, "Delay").text = str(stageInfo[9])
                ET.SubElement(xmlDecimation, "Correction").text = str(stageInfo[9])

            xmlStageGain = ET.SubElement(xmlStage, "StageGain")
            ET.SubElement(xmlStageGain, "Value").text = str(stageInfo[4])
            ET.SubElement(xmlStageGain, "Frequency").text = str(stageInfo[5])

            previousStage = f"{stationInfo[0]}_{stageInfo[0]}_{stageInfo[1]}_{stageInfo[3]}"
            previousXmlStage = xmlStage

    _pretty_print(xmlRoot)

    if len(output) > 0:
        # Save the XML in a file
        byteXml = ET.tostring(xmlRoot)
        with open(output, "wb") as f:
            f.write(byteXml)
        if os.path.isfile(VALIDATOR):
            print(f"Validating {output}")
            subprocess.run(f"java -jar {VALIDATOR} --input {output}", shell=True)
        else:
            print("Validator not found")
    else:
        # Print the XML
        print(ET.tostring(xmlRoot, "unicode"))
