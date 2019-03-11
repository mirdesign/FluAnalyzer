# Minority Variant Calculator
#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# This script is used to calculate minority variants        #
#############################################################

from __future__ import print_function
import locale
import sys
import collections
import os
import glob
import re
from shutil import copyfile

# Load common shared variables and methods
execfile(os.path.join(os.path.dirname(__file__),"_global_constants.py"))
		
		
inputFolder = ""
outputFolder = ""
outputOverViewFolder = ""

# Constants
__MIN_VALUES_PERCENTAGE_CHANGE = 15

__MIN_VAR_BIN_SIZE = 2
__MIN_VAR_BIN_MINIMAL_SIZE = 1
__MIN_VAR_BIN_DIVIDE_HEADER_CHAR = "."
		
__FILE_CONTENT_HEADER_COVERAGE_OK = "CoverageOK"
__FILE_CONTENT_HEADER_VAR_COVERAGE_OK = "VarCovOK"
__FILE_CONTENT_HEADER_COUNT = "Count"
__FILE_CONTENT_HEADER_POSITION = "position"
__FILE_CONTENT_HEADER_DATA = "DATA"

__OUTPUT_OVERVIEW_FILENAME = "overview.csv"
__OUTPUT_MINIMAL_NUMBER_OF_ABOVE_ZERO = 2	
	


# Show some startup information
def ShowSplash():
	Show ("")
	Show ("--------------------------------")
	Show ("Tool: Minority Variant Calculator")
	Show ("Take the output of Variant Detector and calculate the 'Minority variant frequency' per segment per timepoint/location.")
	Show ("version: 0.3")
	Show ("Contact: m.pater@amc.nl | mpater@mirdesign.nl - Maarten Pater")
	Show ("--------------------------------")
	Show ("Usage: python minority_variant_calculator.py folder_with_output_variant_detector folder_output [bin_size=" + str(__MIN_VAR_BIN_SIZE) + "]")
	Show ("Please provide the primer_incl folders from variant_detector.")
	Show ("Expected Filename strategy: [CUSTOM]_[CUSTOM]_[CUSTOM]_[CUSTOM]_[segment]_[CUSTOM].csv")
	Show ("--------------------------------")
		
# Detect header locations
def GetColumnPositionFromHeader(line, SearchForHeader):
	pos = line.find(__INPUT_ROW_SEPERATOR + SearchForHeader)

	# Check if we found the headers
	if (pos == -1):
		QuitAndShowError([	("(GetColumnPositionFromHeader) " + SearchForHeader + " value could not be detected in the headerline."),
							("Provided: " + line),
							("Searched for: " + SearchForHeader)])
							
	pos = pos + len(__INPUT_ROW_SEPERATOR) # Exclude __INPUT_ROW_SEPERATOR char					
	return line[:pos].count(__INPUT_ROW_SEPERATOR)
	
	
# Detect header locations
def GetLastColumnPositionFromHeader(line, SearchForHeader):
	pos = line.find(SearchForHeader)
			
	# Check if we found the headers
	if (pos == -1):
		QuitAndShowError([	("(GetLastColumnPositionFromHeader) " + SearchForHeader + " value could not be detected in the headerline."),
							("Provided: " + line),
							("Searched for: " + SearchForHeader)])
	
	# Now search for the last position
	
	lastPos = pos
	while(not (pos == -1)):
		lastPos = pos
		pos = line.find(SearchForHeader, pos + 1)
	
	return line[:lastPos].count(__INPUT_ROW_SEPERATOR)
	


# Build CSV header
def buildCSVHeaderOverview(file, countHeaders):
	file.write("ID")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("segment")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("position")
	file.write(__OUTPUT_ROW_SEPERATOR)
	
	for header in countHeaders:
		file.write(header)
		file.write(__OUTPUT_ROW_SEPERATOR)

	file.write("isStable")
	file.write("\n")

# Build CSV header
def buildCSVHeader(file):
	file.write("position")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("PropMinVar")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("PropMajVar")
	file.write("\n")

# Save CSV data to file	
# minorityVariantsPerSegment: segment, countHeader, pos, frequency
def saveOverviewCSV(minorityVariantsPerSegment, segmentSizes, countHeaders):
	file = ""
	
	file = open(outputOverViewFolder + "/" + __OUTPUT_OVERVIEW_FILENAME , "w")
	file.write("sep=" + __OUTPUT_ROW_SEPERATOR + "\n")
	buildCSVHeaderOverview(file, countHeaders)

	for keySegment in minorityVariantsPerSegment:
		for pos in range(1, segmentSizes[keySegment] + 1):
			# Check if we have data for this position
			
			fileBuffer = ""
			hasData = False
			numberOfValuesAboveZero = 0
			
			#if(minorityVariantsPerSegment[keySegment][countHeaders[0]].has_key(pos)):
			if(minorityVariantsPerSegment[keySegment][minorityVariantsPerSegment[keySegment].keys()[0]].has_key(pos)):

				fileBuffer += (keySegment + "_" + str(pos))
				fileBuffer += (__OUTPUT_ROW_SEPERATOR)
				fileBuffer += (keySegment)
				fileBuffer += (__OUTPUT_ROW_SEPERATOR)
				fileBuffer += (str(pos))
				fileBuffer += (__OUTPUT_ROW_SEPERATOR)

				minValue = -1
				maxValue = -1

				for header in countHeaders:
					if(header in minorityVariantsPerSegment[keySegment]):
						mvps = minorityVariantsPerSegment[keySegment][header][pos]
						if (mvps > 0):
							numberOfValuesAboveZero += 1
							hasData = True

							if (minValue == -1 or mvps < minValue):
								minValue = mvps
							if (minValue == -1 or mvps > maxValue):
								maxValue = mvps
						#if (mvps == 0):
						#	fileBuffer += "Zero" # Used to translate value '0' to NA in R, otherwise you'll get '10' as '1NA' after translation.
						#else:
						#	fileBuffer += (str(mvps))
						fileBuffer += (str(mvps))
						fileBuffer += (__OUTPUT_ROW_SEPERATOR)
					else:
						fileBuffer += (__OUTPUT_NOT_AVAILABLE_VALUE)
						fileBuffer += (__OUTPUT_ROW_SEPERATOR)
						

				if(hasData):
					isExceeding = False
					if(minValue > 1):
						minMaxDifference = ((maxValue - minValue) / minValue) * 100.0
						isExceeding = abs(minMaxDifference) >= __MIN_VALUES_PERCENTAGE_CHANGE

						if(isExceeding):
							fileBuffer += ("false")
						else:
							fileBuffer += ("true")
							
					elif(minValue >= 0):
						fileBuffer += ("true")
					else:
						fileBuffer += ("false")

						
				fileBuffer += ("\n")
			
			hasEnoughValue = not(minValue == 0 and maxValue == 0)
			
			if(hasData and numberOfValuesAboveZero >= __OUTPUT_MINIMAL_NUMBER_OF_ABOVE_ZERO):
				file.write(fileBuffer)

			

# Save CSV data to file	
# minorityVariantsPerSegment: segment, countHeader, pos, frequency
def saveCSV(minorityVariantsPerSegment, segmentSizes, isBin = False):

	for keySegment in minorityVariantsPerSegment:
		for keyCountHeader in minorityVariantsPerSegment[keySegment]:
			file = ""
			if (isBin):
				file = open(outputFolder + "/" + keySegment + "-bin." + keyCountHeader + ".csv", "w")
			else:
				file = open(outputFolder + "/" + keySegment + "-" + keyCountHeader + ".csv", "w")
			file.write("sep=" + __OUTPUT_ROW_SEPERATOR + "\n")
			
			buildCSVHeader(file)
			
			writePosToFile = {}
			for keyPos in minorityVariantsPerSegment[keySegment][keyCountHeader]:
				writePosToFile[keyPos] = minorityVariantsPerSegment[keySegment][keyCountHeader][keyPos]
						

			for pos in range(1, segmentSizes[keySegment] + 1):
				file.write(str(pos ))
				file.write(__OUTPUT_ROW_SEPERATOR)
			
				if(pos in writePosToFile):
					file.write(str(writePosToFile[pos]))
					file.write(__OUTPUT_ROW_SEPERATOR)
					file.write(str(100.0 - writePosToFile[pos]))
				else:
					file.write(__OUTPUT_NO_DATA)
					file.write(__OUTPUT_ROW_SEPERATOR)
					file.write(__OUTPUT_NO_DATA)
					
				
				file.write("\n")
					
			file.close()



def ReadDataFromFile(file, fileFolder):
	headerNames = []

	skippedFirstLine = False
	headerLine = True
	
	colFirstCountPoint = -1
	colFirstVarCovOKPoint = -1
	colFirstCoverageOKPoint = -1
	# colCoverageOK = -1
	colProteinPosition = -1
	colDATA = -1
	numberOfCountPoints = -1
	segmentSize = -1
	
	# pos, data, count[per CountHeader]
	retVal = {}
	
	for line in open(fileFolder + "/" + file):
		if(not skippedFirstLine): # Skip header
			skippedFirstLine = True
			continue
		
		if(headerLine):		
			
			# Get start positions of counts and the quantity
			colFirstCountPoint = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_COUNT)
			colLastCountPoint =  GetLastColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_COUNT)
			numberOfCountPoints = colLastCountPoint - colFirstCountPoint + 1

			
			# Get start positions of varcoverage
			colFirstVarCovOKPoint = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_VAR_COVERAGE_OK)
			
			# Get start positions of CoverageOK
			colFirstCoverageOKPoint = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_COVERAGE_OK)
			
			colProteinPosition = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_POSITION)
			colDATA = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_DATA)
			
			# get all header names
			for count in range(numberOfCountPoints):
				headerNames.append(GetDataFromFileContent(line, colFirstCountPoint + count))
				
			headerLine = False
			continue
		
		
		# Only accept rows with coverageOk
		cProteinPosition = int(GetDataFromFileContent(line, colProteinPosition))
		segmentSize = cProteinPosition
		
		cData = GetDataFromFileContent(line, colDATA)

		if (not(cProteinPosition in retVal)):
			retVal[cProteinPosition] = {}
		if (not(cData in retVal[cProteinPosition])):
			retVal[cProteinPosition][cData] = {}
	
		for count in range(numberOfCountPoints):
			cVarCovOK = GetDataFromFileContent(line, colFirstVarCovOKPoint + count).lower()
			cCovOK = GetDataFromFileContent(line, colFirstCoverageOKPoint + count).lower()
			
			# Only accept a timepoint/position if its varCovOk and CoverageOk is true
			if (cVarCovOK == "true" and cCovOK == "true"): 
				countValue = GetDataFromFileContent(line, colFirstCountPoint + count).lower()
				if(countValue.isdigit()):
					retVal[cProteinPosition][cData][headerNames[count]] = int(countValue)
	
	return retVal, headerNames, segmentSize
			
			
# rawDataFromFile format: pos, data, count[per CountHeader]
# return: countHeader, pos, minorityVariant 
def GetMinorityVariantPerTimePerPosition (rawDataFromFile, CountHeaders):

	retVal = {}
	for countHeader in CountHeaders: # countHeader = CountT1.. CountT2..
		if (not(countHeader in retVal)):
			retVal[countHeader] = {}
	
		for keyPos in rawDataFromFile:

			counts = []
			for keyData in rawDataFromFile[keyPos]:
				if(countHeader in rawDataFromFile[keyPos][keyData]):	
					counts.append(int(rawDataFromFile[keyPos][keyData][countHeader]))
			
			# Get consensus count
			consensusCount = 0
			if (len(counts) > 0):
				consensusCount = max(counts) 	# Note that if we have 2 data points with equal count (which are the max count in the list), 
												# only 1 will be used as major; due to the fact we only substract max(..) once
			
			# Determine minority variant frequency
			if(consensusCount > 0):
				totalCount = sum(counts)
				frequency = (float(totalCount - consensusCount) / float(totalCount)) * 100.0
				retVal[countHeader][keyPos] = frequency
			else:
				retVal[countHeader][keyPos] = -1
	return retVal	
	
				
# rawDataFromFile format: pos, data, count[per CountHeader]
# return: countHeader-combination, pos, minorityVariant 
def GetMinorityVariantPerBinPerPosition (rawDataFromFile, CountHeaders):

	retVal = {}
	inBin = 0
	binHeader = ""
	
	# Determine bin headers and reserve their present
	for countHeader in CountHeaders:
		inBin = inBin + 1
		if (inBin > __MIN_VAR_BIN_SIZE):
			inBin = 1
			retVal[binHeader] = {} # Store created bin header
			binHeader = ""
		binHeader = binHeader + __MIN_VAR_BIN_DIVIDE_HEADER_CHAR + countHeader
		
	if (inBin >= __MIN_VAR_BIN_MINIMAL_SIZE):
		retVal[binHeader] = {} # Store created bin header if it did not reach the minimal bin size (i.e. for minimal value 2, a list of 5 timepoints will result in a 2, 2, 1 sized bin)
	
		
	for keyPos in rawDataFromFile:
		inBin = 0
		consensusCount = 0
		countsTotal = []
		binHeader = ""
		for countHeader in CountHeaders: # countHeader = CountT1.. CountT2..
			inBin = inBin + 1
			# Build header bin
			if (inBin > __MIN_VAR_BIN_SIZE):
				inBin = 1
				binHeader = ""
				countsTotal = [] # Reset value
				consensusCount = 0 # Reset value
			binHeader = binHeader + __MIN_VAR_BIN_DIVIDE_HEADER_CHAR + countHeader
			
			# Save scores
			counts = []
			for keyData in rawDataFromFile[keyPos]:
				if(countHeader in rawDataFromFile[keyPos][keyData]):	
					counts.append(int(rawDataFromFile[keyPos][keyData][countHeader]))
		
			# Get consensus count of all CountHeaders
			
			# Get consensus count
			maxCount = 0
			if (len(counts) > 0):
				maxCount = max(counts)
				
			consensusCount = consensusCount + maxCount
			countsTotal.extend(counts)
			
			# If we have an existing header-bin-combination, calculate the value with the saved scores
			# Todo: If we had a bin with Size 0, we didn't reserve it. So it won't be found, but we still calculate all values; bit ugly, may resolve this.
			if (binHeader in retVal):
				# Determine minority variant frequency
				if(consensusCount > 0):
					totalCount = sum(countsTotal)
					frequency = (float(totalCount - consensusCount) / float(totalCount)) * 100.0
					retVal[binHeader][keyPos] = frequency
				else:
					retVal[binHeader][keyPos] = -1

	return retVal	
		
# Process found files
# return: segment, countHeader, pos, frequency
def GetVariantPerSegment(Files, fileFolder):
	if(len(Files) == 0):
		Show ("No files to process. Quit program")
		sys.exit(1)
			
	retValuePerTime = {}
	retValuePerBin = {}
	segmentSizes = {}
	retCountHeaders = []
	
	for file in Files:
		segment  = GetFileSegment(file)
		Show("Analysing file: %s" % file)
		
		
		rawDataFromFile, CountHeaders, segmentSize =  ReadDataFromFile(file, fileFolder) # pos, data, count[per time]
		minorityVariantPerTimePerPosition =  GetMinorityVariantPerTimePerPosition(rawDataFromFile, CountHeaders) # countHeader, pos, minorityVariant 
		minorityVariantPerBinPerPosition =  GetMinorityVariantPerBinPerPosition(rawDataFromFile, CountHeaders) # countHeader-bin, pos, minorityVariant 
		
		# Remember the headers
		# todo: Order might be different as expected; if first header is missing in first file, the only existing header will be the first (while this is the second in all other files which do contain this first missing header)
		for header in CountHeaders:
			if not(header in retCountHeaders):
				retCountHeaders.append(header)
		
		retValuePerTime[segment] = minorityVariantPerTimePerPosition
		retValuePerBin[segment] = minorityVariantPerBinPerPosition
		segmentSizes[segment] = segmentSize
		
	return retValuePerTime, retValuePerBin, segmentSizes, retCountHeaders
			
# Write log file
def WriteLogFile(Files):	
		
	file = open(outputFolder + "/" + IdentifiableExperimentName + ".log", "w")
	file.write("Log file for minority variant calculator\n")
	file.write("\n")
	file.write("SETTINGS\n")
	file.write("showMessages: " + ("TRUE" if showMessages else "FALSE")+ "\n")
	file.write("deepDebug: " + ("TRUE" if deepDebug else "FALSE") + "\n")
	file.write("tabularDivide: " + tabularDivide + "\n") 
			
	file.write("inputFolder: " + inputFolder + "\n")
	file.write("outputFolder: " + outputFolder + "\n")
	file.write("outputOverViewFolder: " + outputOverViewFolder + "\n")
	
	file.write("__FILE_VARIANT_OVERVIEW_FILE_EXTENTION: " + str(__FILE_VARIANT_OVERVIEW_FILE_EXTENTION)+ "\n")
	file.write("__OUTPUT_ROW_SEPERATOR: " + str(__OUTPUT_ROW_SEPERATOR)+ "\n")
			
			
	file.write("\n")
	file.write("FILES\n")
	for f in Files:
		file.write(inputFolder + "/" + f + "\n")
		
	file.close()
		
# Initialise arguments
def Initialise():
	# Write numbers with "," instead of "."
	# locale.setlocale(locale.LC_ALL, 'nl_NL')
			
	args = len(sys.argv)
	if(not (args >= 4)):
		Show("ERROR Not all arguments are provided")
		sys.exit(1)
				
	global inputFolder
	global outputFolder
	global outputOverViewFolder
	global __MIN_VAR_BIN_SIZE
			
	inputFolder = sys.argv[1]
	outputFolder = sys.argv[2]
	outputOverViewFolder = sys.argv[3]
	
	if(args > 4):
		__MIN_VAR_BIN_SIZE = sys.argv[4]
		
# Removing any file which isn't AA/nt typed
def RemoveNonSegmentedOverviewFiles(Files):
			
	global DataTypeFound
	remainingFiles = []
				
	for file in Files:
		if((file[-3:] == __FILE_VARIANT_OVERVIEW_FILE_EXTENTION)): # Skip log files, probably created by reading_frame_selector.py
			remainingFiles.append(file)
			
	Show("Found %s typed files: %s" % (__FILE_VARIANT_OVERVIEW_FILE_EXTENTION, str(len(remainingFiles))))
	return remainingFiles
	
		
		
# Main method
def Main():
		
	ShowSplash()
	Initialise()
			
	Show("")
	Show("Searching for variant overview files...")
	Files = GetFileDirectPaths(inputFolder)
	Files = RemoveNonSegmentedOverviewFiles(Files)


	Show("")
	Show ("Calculate minority variant frequency...")
	minorityVariantsPerSegmentPerTime, minorityVariantsPerSegmentPerBin, segmentSizes, CountHeaders = GetVariantPerSegment(Files, inputFolder) # segment, countHeader, pos, frequency
		
		
	if (len(minorityVariantsPerSegmentPerTime) == 0):
		Show("No calculations made since no value met CoverageOk and VarCovOk.")
	else:
		# (key, value)
		# [segment, variants]
		Show("Write output file per time...")
		saveCSV(minorityVariantsPerSegmentPerTime, segmentSizes)
		Show("Write output file per bin...")
		saveCSV(minorityVariantsPerSegmentPerBin, segmentSizes)
		Show("Write overview output...")
		saveOverviewCSV(minorityVariantsPerSegmentPerTime, segmentSizes, CountHeaders)

		
				
	Show("Ready!") 
		
Main()