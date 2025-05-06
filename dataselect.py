#!/opt/antelope/python3.6.5/bin/python

#==================#
# INTRODUCTION


#==================#
# 	TODO
#
# - Write an introduction


#==================#
# PACKAGES

# Generic Packages
import os
import sys
import signal
import warnings
import argparse
import configparser
from typing import List, Tuple
from datetime import datetime, date, timedelta
from dateutil import rrule 
from serviceModule import iso8601DateTime, seedIdentifier
from serviceModule import filterTable
warnings.filterwarnings("ignore")
signal.signal(signal.SIGINT, signal.SIG_DFL)

# Antelope packages
sys.path.append(os.environ['ANTELOPE'] + "/data/python")
from antelope import datascope

#==================#
# DEFINITIONS

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), "config.ini"))
DB_PATH = config.get("DEFAULT", "DB_PATH")
DB_NAME_BASE = config.get("DEFAULT", "DB_NAME_BASE")
DELAY_IN_DAYS = int(config.get("DEFAULT", "DELAY_IN_DAYS"))
PROGRAM_VERSION = config["PROGRAM_INFO"]["DATASELECT_VERSION"]

#==================#
# ARGUMENTS

parser = argparse.ArgumentParser(description='TODO', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Temporal Constraints
startTime = None
endTime = None
parser.add_argument('--starttime', type=iso8601DateTime, help="Starttime of request")
parser.add_argument('--endtime', type=iso8601DateTime, help="Endtime of request")

# SEED Identifiers
networks = []
stations = []
locations = []
channels = []
parser.add_argument('--network', type=seedIdentifier, help="SEED Network code. Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)")
parser.add_argument('--station', type=seedIdentifier, help="SEED Station code. Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)")
parser.add_argument('--location', type=seedIdentifier, help="SEED Location code (use '--' for blank location identifiers). Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)")
parser.add_argument('--channel', type=seedIdentifier, help="SEED Channel code. Multiple items may also be retrieved using a comma separated list. Wildcards can be used ('?' for single characters, '*' for none or multiple characters)")

# Version
parser.add_argument('-V', '--version', help='Show the version of the program and exit', action='version', version=PROGRAM_VERSION)
#parser.add_argument('-V', '--version', help='Show the version of the program and exit', action='version', version='%(prog)s {0}'.format(PROGRAM_VERSION))

# Other
returnPaths = False
parser.add_argument('-P', '--returnpaths', help='Instead of the content, return the paths to the files', action='store_true')

# Parse the arguments
args = parser.parse_args()
#args = parser.parse_args(['--station', 'CAN', '--network', 'IT'])

if(args.starttime):
    startTime = args.starttime
    if (startTime.date() > date.today()-timedelta(days=DELAY_IN_DAYS)):
        sys.exit()
else:
    sys.exit()
# if(args.endtime):
#     endTime = args.endtime
endTime = startTime

dts = []
for dt in rrule.rrule(rrule.DAILY, dtstart=startTime, until=endTime):
    dts.append((str(dt.year)[-2:], dt.month))
Y_Ms = sorted(list(set(dts)))

if(args.network):
    networks = args.network
# networks = [(0, "RF")]
if(args.station):
    stations = args.station[:1]
else:
    sys.exit()
if(args.location):
    locations = []
    # The value '--' is used when the location is blank. So we
    # must remove all entries with that value and leave them blank too
    # This doesn't happen for RegExes
    for location in args.location:
        if(location[0] == False):
            # '--' can't be normally passed directly, as it is otherwise interpreted as a broken argument
            # So we must check any combination of quotes
            if(location[1] == '--' or location[1] == "'--'" or location[1] == '"--"'):
                locations.append((False, ''))
            else:
                # Not a blank location
                locations.append(location)
        else:
            # Is a RegEx
            locations.append(location)
if(args.channel):
    channels = args.channel

if(args.returnpaths):
    returnPaths = args.returnpaths

#==================#
# IMPLEMENTATION

# Read database and retrieve all file paths
files = []
for Y, M in Y_Ms:
    try:
        if not os.path.exists(os.path.join(DB_PATH, f"{DB_NAME_BASE}_{Y}_{M:02d}")):
            continue
        with datascope.closing(datascope.dbopen(os.path.join(DB_PATH, f"{DB_NAME_BASE}_{Y}_{M:02d}"), perm='r')) as db:
            # Prepare the table structure
            # Filter the results as-you-go, to reduce the amount of data
            # between each step
            networkTable = db.lookup(table='network')
            if(len(networks) > 0):
                networkTable = filterTable(networkTable, 'network.net', networks)

            joinAffiliationTable = networkTable.process(['dbjoin affiliation']).subset('affiliation.net == network.net')
            if(len(stations) > 0):
                joinAffiliationTable = filterTable(joinAffiliationTable, 'affiliation.sta', stations)

            joinSiteTable = joinAffiliationTable.process(['dbjoin site']).subset('site.sta == affiliation.sta')
            joinSiteChannelTable = joinSiteTable.process(['dbjoin sitechan']).subset('site.sta == sitechan.sta')
            joinSChanLocTable = joinSiteChannelTable.process(['dbjoin schanloc']).subset('sitechan.chan == schanloc.chan')
            if(len(channels) > 0):
                joinSChanLocTable = filterTable(joinSChanLocTable, 'schanloc.fchan', channels)
            if(len(locations) > 0):
                joinSChanLocTable = filterTable(joinSChanLocTable, 'schanloc.loc', locations)

            table = joinSChanLocTable.process(['dbjoin wfdisc']).subset('wfdisc.chan == schanloc.chan && wfdisc.sta == site.sta')
            if(startTime):
                table = table.subset('wfdisc.jdate >= {0}'.format(startTime.strftime('%Y%j')))
            if(endTime):
                table = table.subset('wfdisc.jdate <= {0}'.format(endTime.strftime('%Y%j')))
            for tablerec in table.iter_record():
                filename = tablerec.filename()
                if(filename[0] == -1 or filename[0] == 1):
                    files.append(filename[1])
    except Exception as exc:
        print(exc)
    
# Write all distinct filepaths
for filename in list(set(files)):
    if(returnPaths):
        print(filename)
    
    # else:
    # 	with open(filename, 'rb') as file:
    # 		data = file.read()
    # 		print(data.hex())
