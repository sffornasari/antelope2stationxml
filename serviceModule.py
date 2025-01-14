#==================#
# INTRODUCTION


#==================#
# 	TODO
#
# - Understand how to coalesce the results. It seems that if all entries have the same URL,
#   then a wildcard '*' can be used. (I.E.: all stations in a netwowk have the same domain
#   for a specific channel. The result for the post would be something like: {net} * * {chan} 20xx-xx-xx 20xx-xx-xx)
# - Write an introduction


#==================#
# PACKAGES

# Generic Packages
import os
import re
import sys
import enum
import json
import argparse
from types import SimpleNamespace
import dateutil.parser
import urllib.parse as urlparse
import xml.etree.ElementTree as ET
from json import JSONEncoder
from typing import Any, List, Tuple
from datetime import datetime, timedelta
from itertools import groupby
from lib2to3.pgen2.parse import ParseError

# Antelope packages
sys.path.append(os.environ['ANTELOPE'] + "/data/python")
from antelope.datascope import *


#==================#
# CLASSES

class RoutingInfo:
    '''
    Contains information used for Routing requests

    Properties:
    - *network*: a ``str`` with the network code
    - *station*: a ``str`` with the station code
    - *location*: a ``str`` with the location code
    - *channel*: a ``str`` with the channel code
    - *domain*: a ``str`` with the valid domain
    - *startDate*: a ``datetime``. The domain is valid from this date onwards
    - *endDate*: a ``datetime``. The domain is valid from this date backwards
    '''

    def __init__(self, network: str, station: str, location: str, channel: str, domain: str, startDate: datetime, endDate: datetime) -> None:
        '''
        Initialize a new ``RoutingInfo``

        *network* is a ``str`` with the network code

        *station* is a ``str`` with the station code

        *location* is a ``str`` with the location code

        *channel* is a ``str`` with the channel code

        *domain* is a ``str`` with the valid domain

        *startDate* is a ``datetime``. The domain is valid from this date onwards

        *endDate* is a ``datetime``. The domain is valid from this date backwards
        '''
        self.network = network
        self.station = station
        self.location = location
        self.channel = channel
        self.domain = domain
        self.startDate = startDate
        self.endDate = endDate

class PostRequest:
    '''
    The data contained in a row of the body of a POST request
    
    Properties:
    - *network*: a ``str`` with the network code
    - *station*: a ``str`` with the station code
    - *location*: a ``str`` with the location code
    - *domain*: a ``str`` with the valid domain
    - *startDate*: a ``datetime``. The domain is valid from this date onwards
    - *endDate*: a ``datetime``. The domain is valid from this date backwards
    '''

    def __init__(self, network: str, station: str, location: str, channel: str, startDate: datetime, endDate: datetime) -> None:
        '''
        Initialize a new ``PostRequest``

        *network* is a ``str`` with the network code

        *station* is a ``str`` with the station code

        *location* is a ``str`` with the location code

        *channel* is a ``str`` with the channel code

        *startDate* is a ``datetime``. The domain is valid from this date onwards

        *endDate* is a ``datetime``. The domain is valid from this date backwards
        '''
        self.network = network
        self.station = station
        self.location = location
        self.channel = channel
        self.startDate = startDate
        self.endDate = endDate

class JsonEncoder(JSONEncoder):
    '''Custom JSON Encoder for objects'''
    def default(self, o: Any) -> Any:
        return o.__dict__

class Level(enum.Enum):
    '''
    Level

    The depth of the information in the StationXML file
    '''
    network = 0
    '''Network (minimum) level'''
    station = 1
    '''Station level'''
    channel = 2
    '''Channel and location level'''
    response = 3
    '''Stage and calibration (maximum) level'''

class Service(enum.Enum):
    '''
    The name of the module to be called
    '''
    station = 'stationxml.py'
    '''stationxml.py, managing StationXML responses'''
    dataselect = 'dataselect.py'
    '''dataselect.py, managing Dataselect (mSEED) responses'''
    routing = 'routing.py'
    '''routing.py, managing Routing responses'''

class PoleZero:
    '''
    A class containing the information necessary for either a Pole or a Zero
    in StationXML PolesZeros stages

    Properties:
    - *number*: an ``int`` with the index of the pole/zero
    - *real*: a ``float`` with the real part of the pole/zero
    - *realError*: a ``float`` with the error on the real part of the pole/zero
    - *imaginary*: a ``float`` with the imaginary part of the pole/zero
    - *imaginaryError*: a ``float`` with the error on the imaginary part of the pole/zero
    '''
    def __init__(self, number: int, real: float, realError: float, imaginary: float, imaginaryError: float) -> None:
        '''
        Initialize a new PoleZero

        *number* is an ``int`` with the index of the pole/zero

        *real* is a ``float`` with the real part of the pole/zero

        *realError* is a ``float`` with the error on the real part of the pole/zero

        *imaginary* is a ``float`` with the imaginary part of the pole/zero

        *imaginaryError* is a ``float`` with the error on the imaginary part of the pole/zero
        '''
        
        self.number = number
        self.real = real
        self.realError = realError
        self.imaginary = imaginary
        self.imaginaryError = imaginaryError

class PolesZeros:
    '''
    Contains the information relative to the Poles and Zeros of a stage in an instrument calibration

    Properties:
    - *normalizationFactor*: a ``float`` with the normalization factor of the stage
    - *poles*: a ``List`` of ``PoleZero`` with data relative to the poles
    - *zeros*: a ``List`` of ``PoleZero`` with data relative to the zeros
    '''

    def __init__(self, normalizationFactor: float,
        poles: List[PoleZero], zeros: List[PoleZero]) -> None:
        '''
        Initialize a new PolesZeros

        *normalizationFactor* is a ``float`` with the normalization factor of the stage

        *poles* is a ``List`` of ``PoleZero`` with data relative to the poles

        *zeros* is a ``List`` of ``PoleZero`` with data relative to the zeros
        '''
        self.normalizationFactor = normalizationFactor
        self.poles = poles
        self.zeros = zeros

class Decimation:
    '''
    Contains the information relative to the Decimation in the stage of an instrument calibration

    Properties:
    - *inputSampleRate*: an ``int`` with the input sample rate of the stage
    - *factor*: an ``int`` with the conversion factor of the stage
    '''
    def __init__(self, inputSampleRate: int, factor: int) -> None:
        '''
        Initialize a new Decimation
        
        *inputSampleRate* is an ``int`` with the input sample rate of the stage

        *factor* is an ``int`` with the conversion factor of the stage
        '''
        self.inputSampleRate = inputSampleRate
        self.factor = factor

class Coefficient:
    '''
    Contains the information relative to the Coefficient in the stage of an instrument calibration

    Properties:
    - *numerator*: a ``float`` with value of the numerator
    - *error*: a ``float`` with the error of the numerator
    '''
    def __init__(self, numerator: float, error: float) -> None:
        '''
        Initialize a new Coefficient
        
        *numerator*: a ``float`` with value of the numerator

        *error*: a ``float`` with the error of the numerator
        '''
        self.numerator = numerator
        self.error = error

class StationInfo:
    '''
    The data contained in a single step of a stage, necessary to build a StationXML
    '''

    def __init__(self, 
        networkCode, networkName, 
        stationCode, stationName, stationLatitude, stationLongitude, stationElevation, stationStartDate, stationEndDate,
        channelCode, channelForeignCode, locationCode, channelStartDate, channelEndDate, channelDepth, channelDescription, calibrationInstrumentName, calibrationSegType, calibrationSampleRate,
        calibrationChannel, calibrationFrequency, calibrationUnits, calibrationTime, calibrationEndTime, calibrationCalibration, calibrationCalper,
        stageChannel, stageTime, stageEndTime, stageStageId, stageGnom, stageInputUnits, stageOutputUnits, stageSampleRate, stageGainCalibration,  stageDecifac, stageIzero,
        polesZeros, decimation, coefficients) -> None:

        self.networkCode = networkCode
        self.networkName = networkName
        self.stationCode = stationCode
        self.stationName = stationName
        self.stationLatitude = stationLatitude
        self.stationLongitude = stationLongitude
        self.stationElevation = stationElevation
        self.stationStartDate = stationStartDate
        self.stationEndDate = stationEndDate
        self.channelCode = channelCode
        self.channelForeignCode = channelForeignCode
        self.locationCode = locationCode
        self.channelStartDate = channelStartDate
        self.channelEndDate = channelEndDate
        self.channelDepth = channelDepth
        self.channelDescription = channelDescription
        self.calibrationInstrumentName = calibrationInstrumentName
        self.calibrationSegType = calibrationSegType
        self.calibrationSampleRate = calibrationSampleRate
        self.calibrationChannel = calibrationChannel
        self.calibrationFrequency = calibrationFrequency
        self.calibrationUnits = calibrationUnits
        self.calibrationTime = calibrationTime
        self.calibrationEndTime = calibrationEndTime
        self.calibrationCalibration = calibrationCalibration
        self.calibrationCalper = calibrationCalper
        self.stageChannel = stageChannel
        self.stageTime = stageTime
        self.stageEndTime = stageEndTime
        self.stageStageId = stageStageId
        self.stageGnom = stageGnom
        self.stageInputUnits = stageInputUnits
        self.stageOutputUnits = stageOutputUnits
        self.stageSampleRate = stageSampleRate
        self.stageGainCalibration = stageGainCalibration
        self.stageDecifac = stageDecifac
        self.stageIzero = stageIzero
        self.polesZeros = polesZeros
        self.decimation = decimation
        self.coefficients = coefficients


#==================#
# ARGUMENT UTILS

def iso8601DateTime(dateTime: str) -> datetime:
	'''
	ISO-8601 DateTime

	Converts a string containing an ISO-8601 DateTime into a ``datetime`` object
	
	*dateTime* is a ``str``, containing an ISO-8601 datetime

	Returns the string parsed in a ``datetime``
	Raises ``argparse.ArgumentTypeError`` if dateTime is not a valid ISO-8601 DateTime
	'''

	try:
		return dateutil.parser.isoparse(dateTime)
	except ParseError:
		raise argparse.ArgumentTypeError("Not a valid ISO-8601 DateTime: {0!r}".format(dateTime))

def seedIdentifier(argument: str) -> List[Tuple[bool, str]]:
	'''
	Converts a string with a SEED Identifier argument into a ``list`` of ``tuple``

	*argument* is a ``str`` containing a list of codes, separated by a comma. The codes
	can include wildcards, namely:
	'*' to indicate any character, repeated any times (even 0)
	'?' to indicate only a single character

	``argument`` is split using the delimiter; then each code, if it contains a wildcard,
	is converted into the equivalent RegEx

	All codes are converted to UpperCase, as it is assumed that the database contains
	only uppercase codes

	Returns a ``list`` of ``tuple``, whose first element is a ``bool`` indicating if the code
	is a RegEx or not, and the second element is the code (or the equivalent RegEx) as a ``str``
	Raises ``argparse.ArgumentTypeError`` for any exception
	'''

	try:
		codes = []
		# All codes in the database are assumed to be UpperCase
		for code in argument.upper().split(','):
			# Check if the code contains any wildcard
			# '*' is '0-N characters'
			# '?' is 'one character'
			if(('*' in code) or ('?' in code)):
				# Convert the wildcards into a true regular expression
				regex = code.replace('?','[A-Z0-9]').replace('*','[A-Z0-9]*')
				codes.append((True, "^{}$".format(regex)))
			else:
				codes.append((False, code))
		return codes
	except:
		raise argparse.ArgumentTypeError("Not a valid SEED Identifier: {0!r}".format(argument))

def latitudeFloat(latitude: str) -> float:
	'''
	Converts a string into a ``float``, checking if it has a valid value for a latitude.

	*latitude* is a ``str`` containig a ``float``

	Returns *latitude* parsed into a ``float``
	Raises ``argparse.ArgumentTypeError`` if the *latitude* is not a ``float``, or if it is not between
	-90.0° and 90.0°
	'''

	try:
		if(-90.0 <= float(latitude) <= 90.0):
			return float(latitude)
		else:
			raise argparse.ArgumentTypeError("Not a valid latitude: {0!r}. Value must be between -90.0° and 90.0°".format(latitude))
	except:
		raise argparse.ArgumentTypeError("Not a valid latitude: {0!r}".format(latitude))

def longitudeFloat(longitude: str) -> float:
	'''
	Converts a string into a ``float``, checking if it has a valid value for a longitude.

	*longitude* is a ``str`` containig a ``float``

	Returns *longitude* parsed into a ``float``
	Raises ``argparse.ArgumentTypeError`` if the *longitude* is not a ``float``, or if it is not between
	-180.0° and 180.0°
	'''

	try:
		if(-180.0 <= float(longitude) <= 180.0):
			return float(longitude)
		else:
			raise argparse.ArgumentTypeError("Not a valid longitude: {0!r}. Value must be between -180.0° and 180.0°".format(longitude))
	except:
		raise argparse.ArgumentTypeError("Not a valid longitude: {0!r}".format(longitude))


#==================#
# DATABASE UTILS

def filterTable(table: Dbptr, columnName: str, tupleCodes: List[Tuple[bool, str]]) -> Dbptr:
	'''
	Filters the records in a table using a list of valid values for a column. 
	Equivalent to a SQL WHERE statement.

	*table* is a ``Dbptr`` to a specific table in the database

	*columnName* is a ``str`` with the full name of a column in table (e.g.: 'network.net')

	*tupleCodes* is a ``list(tuple)``, where each ``tuple`` is composed by a ``bool`` that indicates
	if the other element is a RegEx, and a ``str`` containing the code/RegEx

	Returns a ``Dbptr`` to the table with the filtered records
	'''

	# The codes without a regex are directly stored as validCodes
	# This is the translation of the code:
	# Filter and keep all entries of networs whose first element of the tuple is False (so is not a regex).
	# Take the second element of the typle (so, keep only the code).
	# Store the list of codes in a list.
	validCodes = list(map(lambda tupleCode: tupleCode[1], filter(lambda tupleCode: tupleCode[0] == False, tupleCodes)))
	# Keep the regexes separated in another list
	wildcardCodes = list(map(lambda tupleCode: tupleCode[1], filter(lambda tupleCode: tupleCode[0] == True, tupleCodes)))
	if(len(wildcardCodes) > 0):
		# If there is any regex, then iterate through the entire table to find
		# which entries' codes match the regular expressions
		recordCount = table.query(dbRECORD_COUNT)
		for i in range(recordCount):
			table.record = i
			code = table.getv(columnName)[0]
			for wildcardCode in wildcardCodes:
				if(re.search(wildcardCode, code)):
					validCodes.append(code)
					break

	if(len(validCodes) > 0):
		# This is the translation of the code:
		# Store validNetworks in a dictionary, and then convert it in a list again. This removes any duplicate entry.
		# Convert every code into a string formatted as 'net == "{code}"'.
		# The result of all conversions is an array.
		# Join all the elements of the array, using the string ' || ' as separator.
		# Perform the query on the table.
		# Input: ['RF','IX']
		# Resulting filter: 'net == "RF" || net == "IX"'
		return table.subset(' || '.join(['{0} == "{1}"'.format(columnName, code) for code in list(dict.fromkeys(validCodes))]))
	else:
		# I expect this query to not return anything.
		return table.subset('network.net == ""')


#==================#
# APP UTILS

def getMimeType(format: str) -> str:
    '''
    Get the MimeType for a specific format

    *format* is the name of the MimeType. Managed values are: ``appXml``, ``appJson``, ``xml``,
    ``json``, ``get``, ``post``, ``text``, ``mseed``

    Returns the correct MimeType. Defaults to 'text/plain' for unknown
    or unmanaged types
    '''

    switch = {
        'appXml': 'application/xml',
        'appJson': 'application/json',
        'xml': 'text/xml',
        'json': 'text/plain',
        'get': 'text/plain',
        'post': 'text/plain',
        'text': 'text/plain',
        'mseed': 'application/vnd.fdsn.mseed'
    }
    return switch.get(format.lower(), 'text/plain')

def getStationXmlOutput(level: Level, sortedStationInfos: List[SimpleNamespace]) -> str:
    '''
    Generates a string with an XML from an ordered list of StationXML information

    *level* is the ``Level`` with the depth of the data to return

    *sortedStationInfos* is a ``List`` of ``SimpleNamespace`` with the data to generate a StationXML
    ordered by network code, station code, channel code, stage time and stage id. In reality, it's
    equivalent to and ordered ``List`` of ``StationInfo``

    Returns a ``str`` containing an XML. If sortedStationInfos is empty, returns an empty ``str`` 
    '''

    if(len(sortedStationInfos) == 0):
        return ''
    
    previousNetwork = ''		# {Network code}
    previousStation = ''		# {Station Code}_{Station Start Date}
    previousChannel = ''		# {Station Code}_{Channel+Location Code}_{Channel+Location Start Date}
    previousResponse = ''		# {Station Code}_{Channel+Location Code}_{Calibration Start Epoch}
    previousStage = ''			# {Station Code}_{Channel+Location Code}_{Calibration Start Epoch}_{Stage Id}
    previousXmlNetwork = None
    previousXmlStation = None
    previousXmlChannel = None
    previousXmlResponse = None
    
    xmlRoot = ET.Element('FDSNStationXML')
    xmlRoot.set('xmlns', 'http://www.fdsn.org/xml/station/1')
    xmlRoot.set('schemaVersion', '1.1')
    ET.SubElement(xmlRoot, 'Source').text = 'RF'
    # ET.SubElement(xmlRoot, 'Sender').text = 'PLACEHOLDER'
    ET.SubElement(xmlRoot, 'Module').text = 'stationxml.py'
    # ET.SubElement(xmlRoot, 'ModuleURI').text = 'PLACEHOLDER'
    ET.SubElement(xmlRoot, 'Created').text = str(datetime.now())

    for stationInfo in sortedStationInfos:
        try:
            if(previousNetwork != stationInfo.networkCode):
                xmlNetwork = ET.SubElement(xmlRoot, 'Network')
                xmlNetwork.set('code', stationInfo.networkCode)
                ET.SubElement(xmlNetwork, 'Description').text = stationInfo.networkName
                
                previousNetwork = stationInfo.networkCode
                previousXmlNetwork = xmlNetwork

            if(level.value < Level.station.value):
                continue
        except:
            print('Something happened - Network')
            raise

        try:
            if(previousStation !=  "{0}_{1}".format(stationInfo.stationCode, stationInfo.stationStartDate)):
                xmlStation = ET.SubElement(previousXmlNetwork, 'Station')
                xmlStation.set('code', stationInfo.stationCode)
                xmlStation.set('startDate', str(datetime.strptime(str(stationInfo.stationStartDate), '%Y%j')))
                if(stationInfo.stationEndDate > 0):
                    xmlStation.set('endDate', str(datetime.strptime(str(stationInfo.stationEndDate), '%Y%j')))
                ET.SubElement(xmlStation, 'Latitude').text = str(stationInfo.stationLatitude)
                ET.SubElement(xmlStation, 'Longitude').text = str(stationInfo.stationLongitude)
                ET.SubElement(xmlStation, 'Elevation').text = str(stationInfo.stationElevation * 1000)	# Elevation is in Km
                # ET.SubElement(xmlStation, 'Vault').text = 'PLACEHOLDER'
                ET.SubElement(xmlStation, 'CreationDate').text = str(datetime.strptime(str(stationInfo.stationStartDate), '%Y%j'))
                
                xmlStationSite = ET.SubElement(xmlStation, 'Site')
                ET.SubElement(xmlStationSite, 'Name').text = stationInfo.stationName

                previousStation = "{0}_{1}".format(stationInfo.stationCode, stationInfo.stationStartDate)
                previousXmlStation = xmlStation

            if(level.value < Level.channel.value):
                continue
        except:
            print('Something happened - Station')
            raise
        
        try:
            if(previousChannel != "{0}_{1}_{2}".format(stationInfo.stationCode, stationInfo.channelCode, stationInfo.channelStartDate)):
                xmlChannel = ET.SubElement(previousXmlStation, 'Channel')
                xmlChannel.set('code', stationInfo.channelForeignCode)
                xmlChannel.set('locationCode', stationInfo.locationCode)
                xmlChannel.set('startDate', str(datetime.strptime(str(stationInfo.channelStartDate), '%Y%j')))
                if(stationInfo.channelEndDate > 0):
                    xmlChannel.set('endDate', str(datetime.strptime(str(stationInfo.channelEndDate), '%Y%j')))
                ET.SubElement(xmlChannel, 'Latitude').text = str(stationInfo.stationLatitude)
                ET.SubElement(xmlChannel, 'Longitude').text = str(stationInfo.stationLongitude)
                ET.SubElement(xmlChannel, 'Elevation').text = str(stationInfo.stationElevation * 1000)	# Elevation is in Km
                ET.SubElement(xmlChannel, 'Depth').text = str(stationInfo.channelDepth)
                ET.SubElement(xmlChannel, 'SampleRate').text = str(stationInfo.calibrationSampleRate)

                xmlChannelSensor = ET.SubElement(xmlChannel, 'Sensor')
                ET.SubElement(xmlChannelSensor, 'Description').text = stationInfo.calibrationInstrumentName
                ET.SubElement(xmlChannelSensor, 'Type').text = stationInfo.calibrationSegType

                xmlChannelComment = ET.SubElement(xmlChannel, 'Comment')
                ET.SubElement(xmlChannelComment, 'Value').text = stationInfo.channelDescription

                previousChannel = "{0}_{1}_{2}".format(stationInfo.stationCode, stationInfo.channelCode, stationInfo.channelStartDate)
                previousXmlChannel = xmlChannel
                pass

            if(level.value < Level.response.value):
                continue
        except:
            print('Something happened - Channel')
            raise
        
        try:
            if(previousResponse != "{0}_{1}_{2}".format(stationInfo.stationCode, stationInfo.calibrationChannel, stationInfo.calibrationTime)):
                xmlResponse = ET.SubElement(previousXmlChannel, 'Response')

                xmlResponseInstrumentSensitivity = ET.SubElement(xmlResponse, 'InstrumentSensitivity')
                ET.SubElement(xmlResponseInstrumentSensitivity, 'Value').text = str(1 / float(stationInfo.calibrationCalibration))
                ET.SubElement(xmlResponseInstrumentSensitivity, 'Frequency').text = str(stationInfo.calibrationFrequency)

                xmlResponseInputUnits = ET.SubElement(xmlResponseInstrumentSensitivity, 'InputUnits')
                ET.SubElement(xmlResponseInputUnits, 'Name').text = stationInfo.calibrationUnits

                xmlResponseOutputUnits = ET.SubElement(xmlResponseInstrumentSensitivity, 'OutputUnits')
                ET.SubElement(xmlResponseOutputUnits, 'Name').text = 'Counts'	# TODO Not sure. There is no output units.

                previousResponse = "{0}_{1}_{2}".format(stationInfo.stationCode, stationInfo.calibrationChannel, stationInfo.calibrationTime)
                previousXmlResponse = xmlResponse
        except:
            print('Something happened - Response')
            raise
        
        try:
            if(previousStage != "{0}_{1}_{2}_{3}".format(stationInfo.stationCode, stationInfo.stageChannel, stationInfo.stageTime, stationInfo.stageStageId)):
                xmlStage = ET.SubElement(previousXmlResponse, 'Stage')
                xmlStage.set('number', str(stationInfo.stageStageId))

                xmlStageGain = ET.SubElement(xmlStage, 'StageGain')
                ET.SubElement(xmlStageGain, 'Value').text = str(stationInfo.stageGnom)
                ET.SubElement(xmlStageGain, 'Frequency').text = str(stationInfo.calibrationFrequency)

                previousStage = "{0}_{1}_{2}_{3}".format(stationInfo.stationCode, stationInfo.stageChannel, stationInfo.stageTime, stationInfo.stageStageId)

                if(stationInfo.polesZeros):
                    # Poles and Zeros
                    xmlPolesZeros = ET.SubElement(xmlStage, 'PolesZeros')

                    xmlPolesZerosInputUnits = ET.SubElement(xmlPolesZeros, 'InputUnits')
                    inputName = stationInfo.stageInputUnits
                    if(len(inputName) == 0 or inputName == '-'):
                        inputName = 'Counts'	# Hardcoded as 'counts' if none is present
                    ET.SubElement(xmlPolesZerosInputUnits, 'Name').text = inputName

                    xmlPolesZerosOutputUnits = ET.SubElement(xmlPolesZeros, 'OutputUnits')
                    outputName = stationInfo.stageOutputUnits
                    if(len(outputName) == 0 or outputName == '-'):
                        outputName = 'Counts'	# Hardcoded as 'counts' if none is present
                    ET.SubElement(xmlPolesZerosOutputUnits, 'Name').text = outputName

                    ET.SubElement(xmlPolesZeros, 'NormalizationFactor').text = str(stationInfo.polesZeros.normalizationFactor)
                    ET.SubElement(xmlPolesZeros, 'PzTransferFunctionType').text = 'Laplace (Radians/Second)'	# Hardcoded
                    #ET.SubElement(xmlPolesZeros, 'NormalizationFrequency').text = str(stageInfo[10])	        # TODO Not sure
                    ET.SubElement(xmlPolesZeros, 'NormalizationFrequency').text = str(stationInfo.calibrationCalper)

                    number = 0
                    for pole in stationInfo.polesZeros.poles:
                        xmlPole = ET.SubElement(xmlPolesZeros, 'Pole')
                        xmlPole.set('number', str(number))

                        xmlPoleReal = ET.SubElement(xmlPole, 'Real')
                        xmlPoleReal.text = str(float(pole.real))

                        error = float(pole.realError)
                        if(error < 0):
                            xmlPoleReal.set('plusError', str(0.0))
                            xmlPoleReal.set('minusError', str(error))
                        else:
                            xmlPoleReal.set('minusError', str(0.0))
                            xmlPoleReal.set('plusError', str(error))

                        xmlPoleImaginary = ET.SubElement(xmlPole, 'Imaginary')
                        xmlPoleImaginary.text = str(float(pole.imaginary))

                        error = float(pole.imaginaryError)
                        if(error < 0):
                            xmlPoleImaginary.set('minusError', str(error))
                            xmlPoleImaginary.set('plusError', str(0.0))
                        else:
                            xmlPoleImaginary.set('minusError', str(0.0))
                            xmlPoleImaginary.set('plusError', str(error))

                        number = number + 1

                    number = 0
                    for zero in stationInfo.polesZeros.zeros:
                        xmlZero = ET.SubElement(xmlPolesZeros, 'Zero')
                        xmlZero.set('number', str(number))

                        xmlZeroReal = ET.SubElement(xmlPole, 'Real')
                        xmlZeroReal.text = str(float(zero.real))

                        error = float(zero.realError)
                        if(error < 0):
                            xmlZeroReal.set('plusError', str(0.0))
                            xmlZeroReal.set('minusError', str(error))
                        else:
                            xmlZeroReal.set('minusError', str(0.0))
                            xmlZeroReal.set('plusError', str(error))

                        xmlZeroImaginary = ET.SubElement(xmlPole, 'Imaginary')
                        xmlZeroImaginary.text = str(float(zero.imaginary))

                        error = float(zero.imaginaryError)
                        if(error < 0):
                            xmlZeroImaginary.set('minusError', str(error))
                            xmlZeroImaginary.set('plusError', str(0.0))
                        else:
                            xmlZeroImaginary.set('minusError', str(0.0))
                            xmlZeroImaginary.set('plusError', str(error))

                        number = number + 1
                else:
                    xmlDecimation = ET.SubElement(xmlStage, 'Decimation')

                    ET.SubElement(xmlDecimation, 'InputSampleRate').text = str(stationInfo.decimation.inputSampleRate)
                    ET.SubElement(xmlDecimation, 'Factor').text = str(stationInfo.decimation.factor)
                    ET.SubElement(xmlDecimation, 'Offset').text = str(stationInfo.stageIzero)

                    ET.SubElement(xmlDecimation, 'Delay').text = str(stationInfo.stageGainCalibration)       # TODO Not sure
                    ET.SubElement(xmlDecimation, 'Correction').text = str(stationInfo.stageGainCalibration)  # TODO Not sure

                    xmlCoefficients = ET.SubElement(xmlStage, 'Coefficients')
                    ET.SubElement(xmlCoefficients, 'CfTransferFunctionType').text = 'Digital'	# Hardcoded

                    xmlCoefficientsInputUnits = ET.SubElement(xmlCoefficients, 'InputUnits')
                    inputName = stationInfo.stageInputUnits
                    if(len(inputName) == 0 or inputName == '-'):
                        inputName = 'Counts'	# Hardcoded as 'counts' if none is present
                    ET.SubElement(xmlCoefficientsInputUnits, 'Name').text = inputName

                    xmlCoefficientsOutputUnits = ET.SubElement(xmlCoefficients, 'OutputUnits')
                    outputName = stationInfo.stageOutputUnits
                    if(len(outputName) == 0 or outputName == '-'):
                        outputName = 'Counts'	# Hardcoded as 'counts' if none is present
                    ET.SubElement(xmlCoefficientsOutputUnits, 'Name').text = outputName
                    for coefficient in stationInfo.coefficients:
                        xmlCoefficientsNumerator = ET.SubElement(xmlCoefficients, 'Numerator')
                        xmlCoefficientsNumerator.text = str(float(coefficient.numerator))

                        error = float(coefficient.error)
                        if(error < 0):
                            xmlCoefficientsNumerator.set('minusError', str(error))
                            xmlCoefficientsNumerator.set('plusError', str(0.0))
                        else:
                            xmlCoefficientsNumerator.set('minusError', str(0.0))
                            xmlCoefficientsNumerator.set('plusError', str(error))
        except:
            print('Something happened - Stage')
            raise

    return ET.tostring(xmlRoot, 'unicode')

def getQueryOptions(arguments: List[str]) -> List[Tuple[str, str]]:
    '''
    Create a list of query options using the arguments of a request

    *arguments* is the list of arguments received in a request

    Returns a ``List`` of ``Tuple``, whose first element is the name of the argument, and the
    second is the actual argument, e.g.: ('--network', 'RF')
    '''

    options = []
    # Time Constraints
    #--- Simple-Time
    startTime = arguments.get('starttime')
    endTime = arguments.get('endtime')
    start = arguments.get('start')  # Alternative name to starttime
    end = arguments.get('end')      # Alternative name to endtime
    #--- Window-Time
    startBefore = arguments.get('startbefore')
    startAfter = arguments.get('startafter')
    endBefore = arguments.get('endbefore')
    endAfter = arguments.get('endafter')
    # Channel Constraints
    network = arguments.get('network')
    station = arguments.get('station')
    location = arguments.get('location')
    channel = arguments.get('channel')
    net = arguments.get('net')  # Alternative name to network
    sta = arguments.get('sta')  # Alternative name to station
    loc = arguments.get('loc')  # Alternative name to location
    cha = arguments.get('cha')  # Alternative name to channel
    # Geographic Constraints
    #--- Area-Rectangle
    minLatitude = arguments.get('minlatitude')
    maxLatitude = arguments.get('maxlatitude')
    minLongitude = arguments.get('minlongitude')
    maxLongitude = arguments.get('maxlongitude')
    minLat = arguments.get('minlat')    # Alternative to minlatitude
    maxLat = arguments.get('maxlat')    # Alternative to maxlatitude
    minLon = arguments.get('minlon')    # Alternative to minlongitude
    maxLon = arguments.get('maxlon')    # Alternative to maxlongitude
    #--- Area-Circle
    latitude = arguments.get('latitude')
    longitude = arguments.get('longitude')
    minRadius = arguments.get('minradius')
    maxRadius = arguments.get('maxradius')
    lat = arguments.get('lat')  # Alternative to latitude
    lon = arguments.get('lon')  # Alternative to longitude
    # Output Control
    format = arguments.get('format')
    #nodata = arguments.get('nodata')   # Used internally, no need to retrieve it.
    # Service-Specific Constraints
    #--- Station
    level = arguments.get('level')
    includeRestricted = arguments.get('includerestricted')
    # updatedAfter = arguments.get('updatedafter)'          # Not supported yet
    # matchTimeSeries = arguments.get('matchtimeteries)'    # Not supported yet
    #--- Dataselect
    quality = arguments.get('quality')
    minimumLength = arguments.get('minimumlength')
    longestOnly = arguments.get('longestonly')

    if(startTime):
        options.append(('--starttime', "'{}'".format(startTime)))
    if(endTime):
        options.append(('--endtime', "'{}'".format(endTime)))
    if(start):
        options.append(('--starttime', "'{}'".format(start)))
    if(end):
        options.append(('--endtime', "'{}'".format(end)))
    if(startBefore):
        options.append(('--startbefore', "'{}'".format(startBefore)))
    if(startAfter):
        options.append(('--startafter', "'{}'".format(startAfter)))
    if(endBefore):
        options.append(('--endbefore', "'{}'".format(endBefore)))
    if(endAfter):
        options.append(('--endafter', "'{}'".format(endAfter)))
    if(network):
        options.append(('--network', "'{}'".format(network)))
    if(station):
        options.append(('--station', "'{}'".format(station)))
    if(location):
        options.append(('--location', "'{}'".format(location)))
    if(channel):
        options.append(('--channel', "'{}'".format(channel)))
    if(net):
        options.append(('--network', "'{}'".format(net)))
    if(sta):
        options.append(('--station', "'{}'".format(sta)))
    if(loc):
        options.append(('--location', "'{}'".format(loc)))
    if(cha):
        options.append(('--channel', "'{}'".format(cha)))
    if(minLatitude):
        options.append(('--minlatitude', "'{}'".format(minLatitude)))
    if(maxLatitude):
        options.append(('--maxlatitude', "'{}'".format(maxLatitude)))
    if(minLongitude):
        options.append(('--minlongitude', "'{}'".format(minLongitude)))
    if(maxLongitude):
        options.append(('--maxlongitude', "'{}'".format(maxLongitude)))
    if(minLat):
        options.append(('--minlatitude', "'{}'".format(minLat)))
    if(maxLat):
        options.append(('--maxlatitude', "'{}'".format(maxLat)))
    if(minLon):
        options.append(('--minlongitude', "'{}'".format(minLon)))
    if(maxLon):
        options.append(('--maxlongitude', "'{}'".format(maxLon)))
    if(latitude):
        options.append(('--latitude', "'{}'".format(latitude)))
    if(longitude):
        options.append(('--longitude', "'{}'".format(longitude)))
    if(minRadius):
        options.append(('--minradius', "'{}'".format(minRadius)))
    if(maxRadius):
        options.append(('--maxradius', "'{}'".format(maxRadius)))
    if(lat):
        options.append(('--latitude', "'{}'".format(lat)))
    if(lon):
        options.append(('--longitude', "'{}'".format(lon)))
    if(format):
        options.append(('--format', "'{}'".format(format)))
    # Not used in backend
    # if(nodata):
    #     options.append(('--nodata', "'{}'".format(nodata)))
    if(level):
        options.append(('--level', "'{}'".format(level)))
    if(includeRestricted):
        options.append(('--includerestricted', "'{}'".format(includeRestricted)))
    # Not supported yet
    # if(updatedAfter):
    #     options.append(('--updatedafter', "'{}'".format(updatedAfter)))
    # if(matchTimeSeries):
    #     options.append(('--matchtimeteries', "'{}'".format(matchTimeSeries)))
    if(quality):
        options.append(('--quality', "'{}'".format(quality)))
    if(minimumLength):
        options.append(('--minimumlength', "'{}'".format(minimumLength)))
    if(longestOnly):
        options.append(('--longestonly', "'{}'".format(longestOnly)))
    
    return options


#==================#
# ROUTING UTILS

def getServiceUrl(domain: str, service: str) -> str:
    '''
    Returns the URL to use for a specific service

    *domain* is a ``str`` with the base URL.
    *service* is a ``str`` with the desired service. Valid options are:
    ``dataselect``, ``station``, ``wfcatalog``

    Returns a ``str`` with the URL of the service
    '''

    switch = {
        'dataselect': '/fdsnws/dataselect/1/query',
        'station': '/fdsnws/station/1/query',
        'wfcatalog': '/fdsnws/wfcatalog/1/query'
    }

    if(domain.startswith('https://') or domain.startswith('http://')):
        return '{0}{1}'.format(domain, switch.get(service.lower(), '/fdsnws/dataselect/1/query'))
    else:
        return 'http://{0}{1}'.format(domain, switch.get(service.lower(), '/fdsnws/dataselect/1/query'))

def getOutput(format: str, service: str, groupedRoutingInfos: groupby) -> str:
    '''
    Writes the content of the grouped list in a specific format

    *format* is a ``str`` with the desired format. Valid options are:
    ``xml``, ``json``, ``get``, ``post`` (default is ``xml``)

    *service* is a ``str`` with the desired service. Valid options are:
    ``dataselect``, ``station``, ``wfcatalog``

    *groupedRoutingInfos* is a list of ``RoutingInfo``, ordered and grouped by their domain

    Returns a ``str`` with the formatted output.
    '''

    switch = {
        'xml': getXmlOutput,
        'json': getJsonOutput,
        'get': getGetOutput,
        'post': getPostOutput
    }
    return switch.get(format.lower(), getXmlOutput)(service, groupedRoutingInfos)

def getXmlOutput(service: str, groupedRoutingInfos: groupby) -> str:
    '''
    Generates a string with an XML Document from an ordered and grouped list of routing information

    *service* is a ``str`` with the desired service. Valid options are:
    ``dataselect``, ``station``, ``wfcatalog``

    *groupedRoutingInfos* is a list of ``RoutingInfo``, ordered and grouped by their domain

    Returns a ``str`` containing an XML.
    '''

    # Create an XML object and populate it
    xmlRoot = ET.Element('service')
    for domain, infos in groupedRoutingInfos:
        xmlDatacenter = ET.SubElement(xmlRoot, 'datacenter')
        ET.SubElement(xmlDatacenter, 'name').text = service
        ET.SubElement(xmlDatacenter, 'url').text = getServiceUrl(domain, service)

        for routingInfo in infos:
            xmlParams = ET.SubElement(xmlDatacenter, 'params')
            ET.SubElement(xmlParams, 'net').text = routingInfo.network
            ET.SubElement(xmlParams, 'sta').text = routingInfo.station

            loc = routingInfo.location
            if(len(loc) == 0):
                loc = '*'
            ET.SubElement(xmlParams, 'loc').text = loc

            cha = routingInfo.channel
            if(len(cha) == 0):
                cha = '*'
            ET.SubElement(xmlParams, 'cha').text = cha

            # Dates are stored as YYYYjjj (Julian Date), where
            # YYYY is the year
            # jjj  is the day of the year
            # For example, 2018032 is the 32nd day (1st February) of 2018
            #              2020365 is the 365th day (30th December) of 2020 (a bixestile year)
            # We convert them to a Date
            start = routingInfo.startDate
            if(start < 0):
                start = ''
            else:
                start = datetime.strptime(str(start), '%Y%j')
            ET.SubElement(xmlParams, 'start').text = str(start)

            end = routingInfo.endDate
            if(end < 0):
                end = ''
            else:
                end = datetime.strptime(str(end), '%Y%j')
            ET.SubElement(xmlParams, 'end').text = str(end)

            ET.SubElement(xmlParams, 'priority').text = '1'

    # # Save the XML in a file
    # byteXml = ET.tostring(xmlRoot)
    # with open("routing.xml", "wb") as f:
    #     f.write(byteXml)
    
    return ET.tostring(xmlRoot, 'unicode')

def getJsonOutput(service: str, groupedRoutingInfos: groupby) -> str:
    '''
    Generates a string with a JSON from an ordered list of routing information

    *service* is a ``str`` with the desired service. Valid options are:
    ``dataselect``, ``station``, ``wfcatalog``

    *groupedRoutingInfos* is a list of ``RoutingInfo``, ordered and grouped by their domain

    Returns a ``str`` with a JSON.
    '''

    # Create a list of dictionaries
    routingDictionaries = []
    for domain, infos in groupedRoutingInfos:
        # Create a dictionary. It will be converted into a JSON
        routingDictionary = {
            "name": service,
            "url": getServiceUrl(domain, service),
            "params": []
        }

        for routingInfo in infos:
            net = routingInfo.network
            sta = routingInfo.station
            loc = routingInfo.location
            if(len(loc) == 0):
                loc = '*'

            cha = routingInfo.channel
            if(len(cha) == 0):
                cha = '*'

            # Dates are stored as YYYYjjj (Julian Date), where
            # YYYY is the year
            # jjj  is the day of the year
            # For example, 2018032 is the 32nd day (1st February) of 2018
            #              2020365 is the 365th day (30th December) of 2020 (a bixestile year)
            # We convert them to a Date
            start = routingInfo.startDate
            if(start < 0):
                start = ''
            else:
                start = str(datetime.strptime(str(start), '%Y%j').isoformat())

            end = routingInfo.endDate
            if(end < 0):
                end = ''
            else:
                end = str(datetime.strptime(str(end), '%Y%j').isoformat())

            routingDictionary["params"].append({
                "net": net,
                "sta": sta,
                "loc": loc,
                "cha": cha,
                "start": start,
                "end": end,
            })
        routingDictionaries.append(routingDictionary)

    return json.dumps(routingDictionaries)

def getGetOutput(service: str, groupedRoutingInfos: groupby) -> str:
    '''
    Generates a string with a series of possible GET requests from an ordered and grouped list of routing information

    *service* is a ``str`` with the desired service. Valid options are:
    ``dataselect``, ``station``, ``wfcatalog``

    *groupedRoutingInfos* is a list of ``RoutingInfo``, ordered and grouped by their domain

    Returns a ``str`` with the POST requests, already formatted.
    '''

    output = ''
    for domain, infos in groupedRoutingInfos:
        url = getServiceUrl(domain, service)

        for routingInfo in infos:
            params = {}
            if(len(routingInfo.network) > 0):
                params['net'] = routingInfo.network
            
            if(len(routingInfo.station) > 0):
                params['sta'] = routingInfo.station
            
            if(len(routingInfo.location) > 0):
                params['loc'] = routingInfo.location
            
            if(len(routingInfo.channel) > 0):
                params['cha'] = routingInfo.channel

            if(len(str(routingInfo.startDate)) > 0 and routingInfo.startDate > 0):
                params['start'] = str(datetime.strptime(str(routingInfo.startDate), '%Y%j').isoformat())
            
            if(len(str(routingInfo.endDate)) > 0 and routingInfo.endDate > 0):
                params['end'] = str(datetime.strptime(str(routingInfo.endDate), '%Y%j').isoformat())

            # Compose the URL with the parameters
            url_parts = list(urlparse.urlparse(url))
            query = dict(urlparse.parse_qsl(url_parts[4]))
            query.update(params)
            url_parts[4] = urlparse.urlencode(query)

            output = output + urlparse.urlunparse(url_parts) + '\r\n'
    
    return output

def getPostOutput(service: str, groupedRoutingInfos: groupby) -> str:
    '''
    Generates a string with a series of possible POST requests from an ordered and grouped list of routing information

    *service* is a ``str`` with the desired service. Valid options are:
    ``dataselect``, ``station``, ``wfcatalog``

    *groupedRoutingInfos* is a list of ``RoutingInfo``, ordered and grouped by their domain

    Returns a ``str`` with the POST requests, already formatted.
    '''

    output = ''
    for domain, infos in groupedRoutingInfos:
        output = output + getServiceUrl(domain, service) + '\r\n'

        for routingInfo in infos:
            net = routingInfo.network
            if(len(net) == 0):
                net = '*'
            
            sta = routingInfo.station
            if(len(sta) == 0):
                sta = '*'
            
            loc = routingInfo.location
            if(len(loc) == 0):
                loc = '*'
            
            cha = routingInfo.channel
            if(len(cha) == 0):
                cha = '*'

            start = str(routingInfo.startDate)
            if(len(start) == 0):
                start = '*'
            else:
                start = str(datetime.strptime(str(start), '%Y%j').isoformat())
            
            end = str(routingInfo.endDate)
            if(len(end) == 0 or routingInfo.endDate < 0):
                # If there is no end date, use tomorrow
                tomorrow = datetime.now().date() + timedelta(days=1)
                end = str(tomorrow.isoformat())
            else:
                end = str(datetime.strptime(str(end), '%Y%j').isoformat())

            output = output + '{0} {1} {2} {3} {4} {5}\r\n'.format(net, sta, loc, cha, start, end)
        
        output = output + '\r\n'
    
    return output

def formatBodyOption(option: str) -> str:
    '''
    Format the options received in a POST body, using the full variable name (if present)
    and prepending two hyphens ('--')

    *option* is the name of an argument

    Returns the full name of *option* (if present, otherwise uses *option itself*), with two
    hyphens prepending it (e.g.: 'net' returns '--network')
    '''

    switch = {
        'start': 'starttime',
        'end': 'endtime',
        'net': 'network',
        'sta': 'station',
        'loc': 'location',
        'cha': 'channel',
        'minlat': 'minlatitude',
        'maxlat': 'maxlatitude',
        'minlon': 'minlongitude',
        'maxlon': 'maxlongitude',
        'lat': 'latitude',
        'lon': 'longitude',
    }
    return '--{}'.format(switch.get(option.lower(), option.lower()))


#==================#
# PARSERS

def parseBody(body : str) -> Tuple[List[Tuple[str, str]], List[PostRequest], str, str, Level] :
    '''
    Parse the body of a POST request into a list of arguments, shared among all calls, and
    a list of single requests

    *body* is the content of the the POST request

    Returns a ``Tuple`` with these elements: 
     - a ``List`` of ``Tuple``, whose first element is the name of the argument, and the
       second is the actual argument, e.g.: ('--network', 'RF'). The elements derive from the variables at the
       start of the body (with the structure 'key=value') 
    - a ``List`` of ``PostRequest``, where each of them is a non-variable row  in the body (with the structure
      {net} {sta} {loc} {chan} {start} {end})
    - a ``str`` with the output format (xlm, json, post, or get)
    - a ``str`` with the service (station or dataselect)
    - a ``Level`` with the value of the level (network, station, channel or response). Used only in StationXML
    '''

    rows = body.split('\r\n')
    postRequests = []
    sharedOptions = []
    format = ''
    service = ''
    level = Level.station
    for row in rows:
        rowParts = row.split(' ')
        if(len(rowParts) == 1):
            option = rowParts[0].split('=')
            sharedOptions.append((formatBodyOption(option[0]), "'{}'".format(option[1])))
            if(option[0] == 'format'):
                format = option[1]
            elif(option[0] == 'service'):
                service = option[1]
            elif(option[0] == 'level'):
                level = Level[option[1]]
        elif(len(rowParts) == 6):
            postRequests.append(PostRequest(rowParts[0], rowParts[1], rowParts[2], rowParts[3], rowParts[4], rowParts[5]))
        else:
            raise Exception('Body is not in correct format')
    return (sharedOptions, postRequests, format, service, level)

def parseResponse(format: str, responseLines : List[str]) -> Tuple[bool, str]:
    '''
    Parse a list of string lines in a specific format

    *format* is the name of the format that *responseLines* is expected to represent. Available options are:
    'xml', 'json', 'get', 'post'. Defaults to 'xml'

    *resposneLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is the content of *resposneLines* if parsed successfully,
    otherwise it contains the error message
    '''

    switch = {
        'xml': parseXmlResponse,
        'json': parseJsonResponse,
        'get': parseGetResponse,
        'post': parsePostResponse
    }
    return switch.get(format.lower(), parseXmlResponse)(responseLines)

def parseXmlResponse(responseLines : List[str]) -> Tuple[bool, str]:
    '''
    Try to parse a list of strings as if they were an XML document

    *responseLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is the content of *resposneLines* if parsed successfully,
    otherwise it contains the error message
    '''

    if(responseLines[0] == ''):
        # The response is empty
        return (True, '')
    elif(responseLines[0][0] == '<'):
        # If the first character is a '<', we can assume the rest
        # of the response is an XML
        return (True, ''.join(responseLines))
    else:
        # If there's an error, the backend answers with an explanation 
        # on how to use the command and the cause of the problem. 
        # It appears in the second-to-last row of responseLines
        return (False, responseLines[len(responseLines) - 2])

def parseJsonResponse(responseLines : List[str]) -> Tuple[bool, str]:
    '''
    Try to parse a list of strings as if they were a JSON

    *responseLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is the content of *resposneLines* if parsed successfully,
    otherwise it contains the error message
    '''

    if(responseLines[0] == ''):
        # The response is empty
        return (True, '')
    elif(responseLines[0][0] == '{' or responseLines[0][0] == '['):
        # If the first character is a '{' or '[', we can assume the rest
        # of the response is a JSON
        return (True, ''.join(responseLines))
    else:
        # If there's an error, the backend answers with an explanation 
        # on how to use the command and the cause of the problem. 
        # It appears in the second-to-last row of responseLines
        return (False, responseLines[len(responseLines) - 2])

def parseGetResponse(responseLines : List[str]) -> Tuple[bool, str]:
    '''
    Try to parse a list of strings as if they were a list of GET requests

    *responseLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is the content of *resposneLines* if parsed successfully,
    otherwise it contains the error message.
    '''

    if(responseLines[0] == ''):
        # The response is empty
        return (True, '')
    elif(responseLines[0].startswith('http')):
        # If the response starts with a 'http', we can assume the rest
        # of the response is a GET result.
        # Join them, since we must return a single string
        return (True, '\r\n'.join(responseLines))
    else:
        # If there's an error, the backend answers with an explanation 
        # on how to use the command and the cause of the problem. 
        # It appears in the second-to-last row of responseLines
        return (False, responseLines[len(responseLines) - 2])

def parsePostResponse(responseLines : List[str]) -> Tuple[bool, str]:
    '''
    Try to parse a list of strings as if they were a list of POST requests

    *responseLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is the content of *resposneLines* if parsed successfully,
    otherwise it contains the error message.
    '''

    if(responseLines[0] == ''):
        # The response is empty
        return (True, '')
    elif(responseLines[0].startswith('http')):
        # If the response starts with a 'http', we can assume the rest
        # of the response is a POST result
        # Join them, since we must return a single string
        return (True, '\r\n'.join(responseLines))
    else:
        # If there's an error, the backend answers with an explanation 
        # on how to use the command and the cause of the problem. 
        # It appears in the second-to-last row of responseLines
        return (False, responseLines[len(responseLines) - 2])

def parseMseedResponse(responseLines : List[str]) -> Tuple[bool, str]:
    '''
    Try to parse a list of strings as if they were a MiniSEED file represented as a list of
    hexadecimal values

    *responseLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is the content of *resposneLines* if parsed successfully,
    otherwise it contains the error message.
    ''' 

    if(responseLines[0] == ''):
        # The response is empty
        return (True, '')
    elif(responseLines[0].startswith('usage')):
        # If there's an error, the backend answers with an explanation 
        # on how to use the command and the cause of the problem. 
        # It appears in the second-to-last row of responseLines
        return (False, responseLines[len(responseLines) - 2])
    else:
        # This is most probably a long hex string
        return (True, ''.join(responseLines))

def parseRawResponse(responseLines : List[str]) -> Tuple[bool, List[str]]:
    '''
    Try to parse a list of strings checking if they are a list of raw JSONs

    *responseLines* is a ``List`` of ``str``. Each of them is a row in the terminal
    received as a response from the backend

    Returns a ``Tuple``, whose first element is a ``bool`` indicating the success of the parse
    (``True`` if succeded), while the second one is equal *resposneLines* if successful,
    otherwise it contains the error message.
    ''' 

    lines = []
    for responseLine in responseLines:
        if(responseLine == ''):
            # The line is empty
            continue
        elif(responseLine[0] == '{' and responseLine[len(responseLine) -1] == '}'):
            # If the first character is a '{' and the last is a '}', we can assume the rest
            # of the line is a JSON
            lines.append(responseLine)
        else:
            # If there's an error, the backend answers with an explanation 
            # on how to use the command and the cause of the problem. 
            # It appears in the second-to-last row of responseLines
            print(responseLine)
            return (False, [responseLines[len(responseLines) - 2]])

    return (True, lines)
