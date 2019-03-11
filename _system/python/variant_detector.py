#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
#  Variant detector					                        #
#############################################################
import locale
import sys
import collections
import os
import glob
import re
from shutil import copyfile

# Load common shared variables and methods
execfile(os.path.join(os.path.dirname(__file__),"_global_constants.py"))

# Change these default values if required
Minimal_Value_RowCoverage = 100
__VarCovOK_Minimal_Value_Count = 5
__VarCovOK_Minimal_Value_Percentage = 1 # %
Minimal_Value_Percentage_Change = 15 # %
Minimal_Value_Row_Var_Perc = 1 # % -> Is at least 1 VarPerc on this row AT LEAST this % ?
Minimal_Value_Percentage_Difference = 5 # % -> The value we need to see AT LEAST to mark the position as 'Percentage Difference is OK'
Threshold_Value_Percentage_Difference = 1 #% -> The value of varPerc minimal has to be to take the value in account in the 'Percentage Difference is OK' analysis
Threshold_Value_Max_Reached = 5 # % Is Max value AT LEAST this high?

__RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time = 2
__RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location = 1



# Settings
primerFile = ""
inputFolder = ""
outputFolder = ""
logfile_Name_Prefix = ""
IsTimeInFilesUnique = False
IsLocationInFilesUnique = False
FileNames_UniqueValues = []


# Constants 
DataTypeFound = ""
DataTypeFound_DisplayName = ""
TypeAA_DisplayName = "protein"
TypeNT_DisplayName = "segment"

__FILE_NAME_TIME_DIVIDE = "_" # used to get time from filename
__FILE_NAME_TIME_COL = 2 # used to get time from filename
__FILE_NAME_TIME_PRE_SIZE = 1
__FILE_NAME_LOCATION_COL = 3 # used to get location from filename
__FILE_NAME_LOCATION_PRE_SIZE = 0
__FILE_CONTENT_COL_ROW = 4
__FILE_CONTENT_COL_DEPTH = 3
__FILE_CONTENT_COL_CONSENSUS = 2
__FILE_OUTPUT_APPENDIX_NONPRIMERSONLY = ".NonPrimersOnly"
__FILE_OUTPUT_FOLDER_NONPRIMER = "primer_excl"
__FILE_OUTPUT_FOLDER_WITHPRIMER = "primer_incl"
__FILE_OUTPUT_FOLDER_CONSENSUS = "consensus"
__FILE_OUTPUT_FILENAME_CONSENSUS = "consensus.fasta"
__FILE_OUTPUT_FOLDER_LOG = "log"
__DATA_LETTERS_SELECTED = []
__DATA_LETTERS_AA = ["F", "L", "I", "M", "V", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "S", "G", "*"]
__DATA_LETTERS_NT = ["A", "C", "G", "T"]
__FIRST_COUNT_COL_HEADER_NAME_SELECTED = ""
__FIRST_COUNT_COL_HEADER_NAME_AA = "Total_F"
__FIRST_COUNT_COL_HEADER_NAME_NT = "Total_A"



__READ_DATA_TYPE_TIME = "time"
__READ_DATA_TYPE_DEPTH = "depth"
__OUTPUT_LOW_COVERAGE_VALUE = "LowCov" # LC, R does not support other values then 'NA'.
__OUTPUT_LOW_VAR_COVERAGE_VALUE = "LowVarCov" # LC, R does not support other values then 'NA'.


# Show some startup information
def ShowSplash():
	Show ("")
	Show ("--------------------------------")
	Show ("Tool: Variant detector")
	Show ("Purpose: Parse Multiple influenza segment files and create a timetable to detect variants")
	Show ("version: 1.0")
	Show ("Contact: m.pater@amc.nl | mpater@mirdesign.nl - Maarten Pater")
	Show ("--------------------------------")
	Show ("Usage: python variant_detector.py folder_with_reading_frame_files output_folder primer_positions_file logfile_name_prefix [Minimal_Value_RowCoverage=" + str(Minimal_Value_RowCoverage) + "]")
	Show ("Expected Filename strategy: [CUSTOM]_[CUSTOM]_T[0..n]_[LOCATION]_f[READINGFRAME]-[segment]_[" + __DATATYPE_AMINOACIDS + "/" + __DATATYPE_NUCLEOTIDES + "].csv")
	Show ("Automatic " + __DATATYPE_AMINOACIDS + "/" + __DATATYPE_NUCLEOTIDES + " type detection")
	Show ("Automatic sort based on time or location in file name (T[0..n] or NK or TA)")
	Show ("--------------------------------")

	
# Get AA or NA type information
def GetFileDNAProteinType(file):
	return file[file.rfind("_") + 1: file.rfind("_") + 2 + 1]
	

# Removing any file which isn't AA/nt typed
def RemoveNonReadingFramedFiles(Files):
	
	global DataTypeFound
	global DataTypeFound_DisplayName
	remainingFiles = []
	
	bAAFound = False
	bNTFound = False
	
	for file in Files:
		if (GetFileDNAProteinType(file) == __DATATYPE_AMINOACIDS or GetFileDNAProteinType(file) == __DATATYPE_NUCLEOTIDES):
			if (GetFileDNAProteinType(file) == __DATATYPE_AMINOACIDS):
				if(bNTFound): # Do not accept mix of AA and nt typed files
					# ERROR
					QuitAndShowError([("Both AA and nt typed files in input folder found."),
										("Cannot decide which files to use.")])
				bAAFound = True
				DataTypeFound = __DATATYPE_AMINOACIDS
				DataTypeFound_DisplayName = TypeAA_DisplayName
			else:
				if(bAAFound): # Do not accept mix of AA and nt typed files\
					# ERROR
					QuitAndShowError([("Both AA and nt typed files in input folder found."),
									("Cannot decide which file to use.")])
					sys.exit(1)
				bNTFound = True
				DataTypeFound = __DATATYPE_NUCLEOTIDES
				DataTypeFound_DisplayName = TypeNT_DisplayName
			if(not (file[-3:] == "log")): # Skip log files, probably created by reading_frame_selector.py
				remainingFiles.append(file)
	
	Show("Found %s typed files: %s" % (DataTypeFound, str(len(remainingFiles))))
	return remainingFiles

	
# Parse a line of tabular diveded data and gather the AA/nt data
# Return an array of counts in the same order of __DATA_LETTERS_SELECTED
def GetAACountsFromFileLine(line, firstAAPosition):
	aaCounts = []
	for cnt in range (0, len(__DATA_LETTERS_SELECTED)):
		
		content = GetDataFromFileContent(line, firstAAPosition + cnt, tabularDivide)

		if (not (content == "")):
			aaCounts.append(int(content))
		else:
			aaCounts.append(0)
			
	return aaCounts
	
# Parse a line of tabular diveded data and gather the depths per position
def GetDepthFromFileLine(line):
	depth = -1

	content = GetDataFromFileContent(line, __FILE_CONTENT_COL_DEPTH, tabularDivide)

	if (not (content == "")):
		depth = int(content)
	else:
		# ERROR
		QuitAndShowError([("(GetDepthFromFileLine) Depth could not be read from this line. Depth data is empty."),
							("Provided: " + line)])
	
	return depth
	
# Read tabular divided files and retrieve their AA counts
# Return per row, all countings of all files
# Tree:
# pos 1 
#      \-- file 1
#				 \---F:1
#					-L:322
#					-...etc..
#      \-- file 2
#				 \---F:43
#					-L:632
#					-...etc..
# pos 2
#      \-- file 1
#				 \---F:0
#					-L:3
#					-...etc..
#      \-- file 2
#				 \---F:0
#					-L:11
#					-...etc..
#

# Check if all are same Time or same locations


# Get time from file name
def GetUniqueSortingValueFromFileName(fileName, UniqueValueLocation = -1, UniqueValueSize = -1, isInt = True):

	# Check if we need to use the automaticaly detected values
	if (UniqueValueLocation == -1):
		if(IsTimeInFilesUnique):
			UniqueValueLocation = __FILE_NAME_LOCATION_COL
			isInt = False
		else:
			UniqueValueLocation = __FILE_NAME_TIME_COL
			
	# Check if we need to use the automaticaly detected values
	if (UniqueValueSize == -1):
		if(IsTimeInFilesUnique):
			UniqueValueSize = __FILE_NAME_LOCATION_PRE_SIZE
			isInt = False
		else:
			UniqueValueSize = __FILE_NAME_TIME_PRE_SIZE
		
	# File name example: H3N2_V1298_T14_TA_minQ25__f3.2013-12-21-A.Stockholm.40.2013-EPI_ISL_155657-direct-A.H3N2-HA_AA.txt
	
	# Get time from file name
	uniqueValueInFile = GetDataFromFileContent(fileName, UniqueValueLocation, __FILE_NAME_TIME_DIVIDE)
	
	# Remove leading "T" or other values
	uniqueValueInFile = uniqueValueInFile [UniqueValueSize:]
	
	try:
		if (isInt):
			return int(uniqueValueInFile)
		else:
			return uniqueValueInFile
	except ValueError:
		# ERROR
		QuitAndShowError([	("(GetUniqueSortingValueFromFileName) Could not determine TIME from filename."),
							("Provided fileName: " + fileName),
							("Provided UniqueValueLocation: " + str(UniqueValueLocation)),
							("Provided UniqueValueSize: " + str(UniqueValueSize)),
							("Supported: [CUSTOM]_[CUSTOM]_T[0..n]_[__CUSTOM__]-[segment]_[" + __DATATYPE_AMINOACIDS + "/" + __DATATYPE_NUCLEOTIDES + "].ext"),
							("Did you pay good attention to the underscores \"_\" and the \"T\" prefix for time?")])
		
		
# Search for the first AA position in the header line
def GetFirstAAPositionFromHeader(line):
	firstAAPosition = line.find(__FIRST_COUNT_COL_HEADER_NAME_SELECTED)
	
	if (firstAAPosition == -1):
		QuitAndShowError([	("(GetFirstAAPositionFromHeader) Could not find the first AA position based on header information."),
							("Provided: " + line),
							("Searched for: " + __FIRST_COUNT_COL_HEADER_NAME_SELECTED)])
		
	return line[:firstAAPosition].count(tabularDivide)
	
	


# Get all data from the input files
# readDataType: Can be "time" or "depth"
def readDataFromFile(file, similarFiles, readDataType):

	files = []
	files.append(file)
	files.extend(similarFiles)
	
	# Sort, so timelines are correct (assume time is 1st number in file name)
	files = sorted(files, key=lambda filename: [GetUniqueSortingValueFromFileName(filename)])
	
	timeTable = []
	
	for cf in files:
		skippedFirstLine = False
		firstAAPosition = -1
		
		for line in open(inputFolder + "/" + cf):
			if(not skippedFirstLine):
				firstAAPosition = GetFirstAAPositionFromHeader(line)
				skippedFirstLine = True
				continue
			
			cRow = int(GetDataFromFileContent(line, __FILE_CONTENT_COL_ROW, tabularDivide)) - 1 #-1 to let it align with 0 based array calling		
			if  (cRow < 0):
				continue
						
			timeTableLen = len(timeTable)
			if(cRow >= timeTableLen):
				for i in range (timeTableLen, cRow + 1):
					timeTable.append({})

			data = None
			
			if (readDataType == __READ_DATA_TYPE_TIME):
				data = GetAACountsFromFileLine(line, firstAAPosition)
			elif (readDataType == __READ_DATA_TYPE_DEPTH):
				data = int(GetDepthFromFileLine(line))
			else:
				QuitAndShowError([("(readDataFromFile) Reading data type not supported."),
									("Provided: " + readDataType),
									("Supported types: [" + __READ_DATA_TYPE_TIME + "/" + __READ_DATA_TYPE_DEPTH + "]")])
								
			timeTable[cRow][cf] = data

			
	# Show debug info	
	#for i in range (0, len(timeTable)):
	#	for j in range (0, len(timeTable[i])):
	#		show = "[%d][%s] = " % (i,timeTable[i].keys()[j])
	#		for cnt in range (0, len(__DATA_LETTERS_SELECTED)):
	#			show = show + ("(%s:%d) " % (__DATA_LETTERS_SELECTED[cnt], timeTable[i].values()[j][cnt]))
	#DeepDebug(show)
		
	for tb in range(0, len(timeTable)):	
		timeTable[tb] = collections.OrderedDict(sorted(timeTable[tb].items(), key=lambda filename: [GetUniqueSortingValueFromFileName(filename[0])]))
		
	return timeTable

# Has current AA/nt enough depth on current position for current time?
def IsCoverageOK(positionTimeTable):
	timesTotal = 0
	
	# Calculate depth
	for dataI in range(0, len(__DATA_LETTERS_SELECTED)):
		timesTotal = timesTotal + positionTimeTable[dataI]
	
	if (timesTotal >= Minimal_Value_RowCoverage):
		return True
		
	return False
		
# Has current AA/nt enough depth on current position for all time?
def IsRowCoverageOK(timeTable, pos, data):
	
	if(len(timeTable) == 0):
		return False
	if(len(timeTable[pos]) == 0):
		return False
		
	for time in range (0, len(timeTable[pos])):
		if (not(IsCoverageOK(timeTable[pos].values()[time]))):
			return False
		
	return True
	
def IsVarCovOK(timeTable, pos, data, time):
	timesTotalBasePercentage = 0
	cTime = timeTable[pos].values()[time][data]
	
	# Count above threshold?
	if (cTime == 0):
		return True
	if (cTime < __VarCovOK_Minimal_Value_Count):
		return False
		
	# Percentage of current AA/nt seen over total AA/nt in current time (file) on this position
	timesTotal = 0
	cPercentage = 0
	for dataI in range(0, len(__DATA_LETTERS_SELECTED)):
		timesTotal = timesTotal + timeTable[pos].values()[time][dataI]
	if (timesTotal > 0):
		timesTotalBasePercentage = 100.0 / timesTotal
		cPercentage = timesTotalBasePercentage * cTime
	
	# Count related to coverage above threshold?
	if (cPercentage < __VarCovOK_Minimal_Value_Percentage):
		return False
		
	return True
	
# Has current AA/nt always >5 and >1% (based on coverage of current time)
def IsRowMultiVarCovOK(timeTable, pos, data, minimal_value):
	timesTotalBasePercentage = 0
	
	NumberOfOkTimes = 0
	
	for time in range (0, len(timeTable[pos])):
		cCount = timeTable[pos].values()[time][data]
		if (IsVarCovOK(timeTable, pos, data, time) and IsCoverageOK(timeTable[pos].values()[time]) and cCount > 0):
			NumberOfOkTimes = NumberOfOkTimes + 1

	return NumberOfOkTimes >= minimal_value
	
# Has current AA/nt, on current position, always value 0 in all times (times=files)
def IsCurrentAANul(timeTable, pos, aa):
	for time in range (0, len(timeTable[pos])):
		if(timeTable[pos].values()[time][aa] > 0):
			return False	
	return True
	
# Is there just 1 AA/nt on this position with value (so, is this AA/nt stable?)
def IsCurrentTimeSingleAA(timeTable, pos):
	isSingle = False
	lastFoundData = -1
	for time in range (0, len(timeTable[pos])):
		for cData in range(0, len(__DATA_LETTERS_SELECTED)):
			if(timeTable[pos].values()[time][cData] > 0):
				if(isSingle == True and not (lastFoundData == cData)):
					return False
				else:
					lastFoundData = cData
					isSingle = True
	
	return isSingle
	
# Build CSV header
def buildCSVHeader(timeTable, LowestTimePointSeen, HighestTimePointSeen):
	retVal = []
	
	retVal.append("ID")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("segment")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("position")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("DATA")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	# retVal.append("PositionIsConserved")
	# retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("RowCountsZero")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("RowAllCoverageOK")
	retVal.append(__OUTPUT_ROW_SEPERATOR)

	#for i in range(0, dataPoints):
	for i in range(LowestTimePointSeen, HighestTimePointSeen + 1):
		if(not(IsLocationInFilesUnique)):
			# Convert time point to location name"
			retVal.append("CoverageOK" + str(FileNames_UniqueValues[i]))
		else:
			realNumber = str(i)
			zeroSpacings = '0' * int(len(str(HighestTimePointSeen)) - len(realNumber))
			spacedNumber = zeroSpacings + realNumber
			retVal.append("CoverageOKT" + spacedNumber)
		retVal.append(__OUTPUT_ROW_SEPERATOR)
	
	#for i in range(0, dataPoints):
	for i in range(LowestTimePointSeen, HighestTimePointSeen + 1):
		if(not(IsLocationInFilesUnique)):
			# Convert time point to location name"
			retVal.append("VarCovOK" + str(FileNames_UniqueValues[i]))
		else:
			realNumber = str(i)
			zeroSpacings = '0' * int(len(str(HighestTimePointSeen)) - len(realNumber))
			spacedNumber = zeroSpacings + realNumber
			retVal.append("VarCovOK" + spacedNumber)
		retVal.append(__OUTPUT_ROW_SEPERATOR)
		
		
	retVal.append("RowMultiVarCovOK") #Time
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	#retVal.append("RowMultiVarCovOKLocation")
	#retVal.append(__OUTPUT_ROW_SEPERATOR)
	
	retVal.append("IsPrimer")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	
	
	#for i in range(0, dataPoints):
	for i in range(LowestTimePointSeen, HighestTimePointSeen + 1):
		if(not(IsLocationInFilesUnique)):
			# Convert time point to location name"
			retVal.append("Coverage" + str(FileNames_UniqueValues[i]))
		else:
			realNumber = str(i)
			zeroSpacings = '0' * int(len(str(HighestTimePointSeen)) - len(realNumber))
			spacedNumber = zeroSpacings + realNumber
			retVal.append("CoverageT" + spacedNumber)
		retVal.append(__OUTPUT_ROW_SEPERATOR)
		
	#for i in range(0, dataPoints):
	for i in range(LowestTimePointSeen, HighestTimePointSeen + 1):
		if(not(IsLocationInFilesUnique)):
			# Convert time point to location name
			retVal.append("Count" + str(FileNames_UniqueValues[i]))
		else:
			realNumber = str(i)
			zeroSpacings = '0' * int(len(str(HighestTimePointSeen)) - len(realNumber))
			spacedNumber = zeroSpacings + realNumber
			retVal.append("CountT" + spacedNumber)
		retVal.append(__OUTPUT_ROW_SEPERATOR)
		
	#for i in range(0, dataPoints):
	for i in range(LowestTimePointSeen, HighestTimePointSeen + 1):
		if(not(IsLocationInFilesUnique)):
			# Convert time point to location name
			retVal.append("VarPerc" + str(FileNames_UniqueValues[i]))
		else:
			realNumber = str(i)
			zeroSpacings = '0' * int(len(str(HighestTimePointSeen)) - len(realNumber))
			spacedNumber = zeroSpacings + realNumber
			retVal.append("VarPercT" + spacedNumber)
		retVal.append(__OUTPUT_ROW_SEPERATOR)
		
	# Do not show percentage change between 2 time points, since we now use min/max as filter values
	# for i in range(1, numTimes):
		# retVal.append("VarPercChangef" + str(i ) + "f" + str(i + 1))
		# retVal.append(__OUTPUT_ROW_SEPERATOR)
		
		
		
	retVal.append("PotentialMaj")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("PercChngOK")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("PercDiffOk")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("PercDiffLowRegion")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("MaxAtLeast")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("MinMaxFirstVarf")
	retVal.append(__OUTPUT_ROW_SEPERATOR) 
	retVal.append("MinMaxSecondVarf")
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("FirstLastChangeVarf") # was: MinMaxChangeFirstVarf
	retVal.append(__OUTPUT_ROW_SEPERATOR)
	retVal.append("MinMaxChangeVarf") # was: MinMaxChangeSecondVarf
	#retVal.append(__OUTPUT_ROW_SEPERATOR)


	
	return retVal
	
def GetCombinedRowID(Segment, pos, currentAA):
	retVal = ""
	retVal = Segment
	retVal += "_"
	retVal += str(pos + 1)
	retVal += "_"
	retVal += __DATA_LETTERS_SELECTED[currentAA]
	return retVal

# check if current position is an primer position
def IsPrimerPosition(PrimerLocations, Segment, pos):
	
	if (Segment in PrimerLocations):
		if(pos in PrimerLocations[Segment]):
			return True
			
	return False
	

# Input:
# [580][0.HA_AA.txt] = (F:0) (L:0) (I:0) (M:0) (V:0) (P:0) (T:0) (A:0) (Y:0) (H:0) (Q:0) (N:0) (K:231) (D:0) (E:0) (C:0) (W:0) (R:0) (S:0) (G:0) (*:0)
# [580][1.HA_AA.txt] = (F:0) (L:0) (I:0) (M:0) (V:0) (P:0) (T:0) (A:0) (Y:0) (H:0) (Q:0) (N:0) (K:11) (D:0) (E:0) (C:0) (W:0) (R:0) (S:0) (G:0) (*:0)
# [581][0.HA_AA.txt] = (F:0) (L:0) (I:0) (M:0) (V:0) (P:0) (T:181) (A:0) (Y:0) (H:0) (Q:0) (N:0) (K:0) (D:0) (E:0) (C:0) (W:0) (R:0) (S:0) (G:0) (*:0)
# [581][1.HA_AA.txt] = (F:0) (L:0) (I:0) (M:0) (V:0) (P:2) (T:8) (A:0) (Y:0) (H:0) (Q:0) (N:0) (K:0) (D:0) (E:0) (C:0) (W:0) (R:0) (S:0) (G:0) (*:0)
#
# Output:
#				isConsevered	t1	t2	
#	ha_140_d	0				100	95	
#	ha_140_n	0				0	5	
#	ha_141_d	1				100	100	<-- only AA on this position with value > 0
#
# Transform timetable from file to count output	
def transformTimeTableToCSV(timeTable, numTimes, readingDepth, numFiles, LowestTimePointSeen, HighestTimePointSeen, Segment, PrimerLocations, TimePoints, NonPrimersOnly = False):
	retVal = []

	timePoints = -1
	if(not(IsLocationInFilesUnique)):
		timePoints = dataPoints
	else:
		timePoints = numFiles

	retVal.append(buildCSVHeader(timeTable, LowestTimePointSeen, HighestTimePointSeen))


	for pos in range (0, len(timeTable)): # each position

		availableTimes = []
		if(IsLocationInFilesUnique):
			for key in timeTable[pos].keys():
				availableTimes.append(GetUniqueSortingValueFromFileName(key))
		else:
			for key in range(0, len(timeTable[pos]) + 1):
				availableTimes.append(key)
		
		for aa in range (0, len(__DATA_LETTERS_SELECTED)): # each AA
			
			sRowName = GetCombinedRowID(Segment, pos, aa) # Row Name
			isRowCountsZero = IsCurrentAANul(timeTable, pos, aa) # is this AA only 0 over all times?
			isRowCoverageOK = IsRowCoverageOK(timeTable, pos, aa) # is this pos coverage >x over all times?
			isRowMultiVarCovOK_Time = IsRowMultiVarCovOK(timeTable, pos, aa, __RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time) # is this AA >5 and >1% over all times?
			isRowMultiVarCovOK_Location = IsRowMultiVarCovOK(timeTable, pos, aa, __RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location) # is this AA >5 and >1% over all times?
			isPrimerPosition = IsPrimerPosition(PrimerLocations, Segment, pos + 1)
			
		
		
			cRow = []
			# Build row value
			cRow.append(sRowName)
			cRow.append(__OUTPUT_ROW_SEPERATOR)
			cRow.append(Segment)
			cRow.append(__OUTPUT_ROW_SEPERATOR)
			cRow.append(str(pos + 1))
			cRow.append(__OUTPUT_ROW_SEPERATOR)
			cRow.append(__DATA_LETTERS_SELECTED[aa])
			cRow.append(__OUTPUT_ROW_SEPERATOR)
			
			# RowCountsZero
			if(isRowCountsZero):
				cRow.append("true")
			else:
				cRow.append("false")
			cRow.append(__OUTPUT_ROW_SEPERATOR)
			
			# RowCoverageOK
			if(isRowCoverageOK):
				cRow.append("true")
			else: 
				cRow.append("false")
			cRow.append(__OUTPUT_ROW_SEPERATOR)
			
			# Add CoverageOKTx values
			previousTimePoint = 0
			timeCorrection = 0
			fileTime = 0

			
			for time in range (LowestTimePointSeen, HighestTimePointSeen + 1):
				# Write data
				if (not(time <= max(TimePoints))):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
				elif(not(time in availableTimes)):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
					timeCorrection = timeCorrection + 1
				else:
					if(IsCoverageOK(timeTable[pos].values()[fileTime])):
						cRow.append("true")
					else:
						cRow.append("false")
					fileTime = fileTime + 1
				cRow.append(__OUTPUT_ROW_SEPERATOR)				
			
			
		
			# Add VarCovOKTx values
			previousTimePoint = 0
			timeCorrection = 0
			fileTime = 0
			for time in range (LowestTimePointSeen, HighestTimePointSeen + 1):
				# Write data
				if (not(time <= max(TimePoints))):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
				elif(not(time in availableTimes)):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
					timeCorrection = timeCorrection + 1
				else:
					if(IsVarCovOK(timeTable, pos, aa, fileTime)):
						cRow.append("true")
					else:
						cRow.append("false")
					fileTime = fileTime + 1
				cRow.append(__OUTPUT_ROW_SEPERATOR)		
				
				
			if(IsLocationInFilesUnique):
				# RowMultiVarCovOK Time (RowMultiVarCovOKTime)
				if(isRowMultiVarCovOK_Time):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
			else:
				# RowMultiVarCovOK Location (RowMultiVarCovOKTime)
				if(isRowMultiVarCovOK_Location):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
					
			# IsPrimer	
			if(isPrimerPosition):
				cRow.append("true")
			else:
				cRow.append("false")
			cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				
			# Add coverageTx
			previousTimePoint = 0
			timeCorrection = 0
			fileTime = 0
			for time in range (LowestTimePointSeen, HighestTimePointSeen + 1):
				# Write data
				if (not(time <= max(TimePoints))):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
				elif(not(time in availableTimes)):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
					timeCorrection = timeCorrection + 1
				else:
					cRow.append(readingDepth[pos].values()[fileTime])
					fileTime = fileTime + 1
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				
			# Add countTx
			previousTimePoint = 0
			timeCorrection = 0
			fileTime = 0

			for time in range (LowestTimePointSeen, HighestTimePointSeen + 1):
				# Write data
				if (not(time <= max(TimePoints))):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
				elif(not(time in availableTimes)):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
					timeCorrection = timeCorrection + 1
				else:
					cRow.append(timeTable[pos].values()[fileTime][aa])
					fileTime = fileTime + 1
				cRow.append(__OUTPUT_ROW_SEPERATOR)
	
	

			# Add varPercTx %
			varPercMinInit = 10000
			varPercMaxInit = -10000
			varPercMin = varPercMinInit
			varPercMax = varPercMaxInit
			
			varPecMin_DefaultValue = varPercMin
			varPercMax_DefaultValue = varPercMax
			varPercMinSetLast = False
			
			
			previousTimePoint = 0
			timeCorrection = 0
			fileTime = 0
			for time in range (LowestTimePointSeen, HighestTimePointSeen + 1):
				# Write data
				if (not(time <= max(TimePoints))):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
				elif(not(time in availableTimes)):
					cRow.append(__OUTPUT_NOT_AVAILABLE_VALUE)
					timeCorrection = timeCorrection + 1
				else:
					cVarCovOK = IsVarCovOK(timeTable, pos, aa, fileTime)
					cCoverageOk = IsCoverageOK(timeTable[pos].values()[fileTime])
					if (not(cCoverageOk)):
						cRow.append(__OUTPUT_LOW_COVERAGE_VALUE)
					elif (not(cVarCovOK)):
						cRow.append(__OUTPUT_LOW_VAR_COVERAGE_VALUE)
					else:
						value = float(timeTable[pos].values()[fileTime][aa])
						depth = float(readingDepth[pos].values()[fileTime])
						
						if (value == 0 or depth == 0):
							cRow.append(0)
							perc = 0
						else:
							perc = (value * 100.0) / depth
							cRow.append(str(locale.format('%f', perc)))
							
						if (perc < varPercMin):
							varPercMin = perc
							varPercMinSetLast = True
							
						if (perc > varPercMax):
							varPercMax = perc
							varPercMinSetLast = False
					fileTime = fileTime + 1
							
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
			# Add time difference %
			lengthDataAvailable = len(timeTable[pos])
		
			IsChangeExceeding = False
			IsPercDifferenceExceeding = False
			varPercChangeFirstLast_First_Init = -1
			varPercChangeFirstLast_Last_Init = -1
			varPercChangeFirstLast_First = varPercChangeFirstLast_First_Init
			varPercChangeFirstLast_Last = varPercChangeFirstLast_Last_Init
			for time in range (LowestTimePointSeen, HighestTimePointSeen + 1):
				if (not(time >= lengthDataAvailable)):
					# Remember our first percentage value
					if (varPercChangeFirstLast_First == -1):
						if (IsVarCovOK(timeTable, pos, aa, time - 1)):
							firstDepth = float(readingDepth[pos].values()[time - 1])
							firstPercentage = 0.0
							firstData = float(timeTable[pos].values()[time - 1][aa])

							if (not(firstData == 0 or firstDepth == 0)):
								firstPercentage = (firstData * 100.0) / firstDepth
							varPercChangeFirstLast_First = firstPercentage
						
					
					if (IsVarCovOK(timeTable, pos, aa, time)):
						secondDepth = float(readingDepth[pos].values()[time])
						secondPercentage = 0.0
						secondData = float(timeTable[pos].values()[time][aa])
						if (not(secondData == 0 or secondDepth == 0)):
							secondPercentage = (secondData * 100.0) / secondDepth
						
						# Remember our last percentage value
						if (time == timePoints - 1):
							varPercChangeFirstLast_Last = secondPercentage
					else:
						# Last is LowCount, so reset to init
						varPercChangeFirstLast_Last = varPercChangeFirstLast_Last_Init 
					
							
		
			if(varPercMin == varPercMinInit and varPercMax == varPercMaxInit):
				cRow.append("false") #PotentialMaj
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append("false") #PercChngOK
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append("false") #PercDiffOk
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append("false") #PercDiffLowRegion
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append("false") # MaxAtLeast
				cRow.append(__OUTPUT_ROW_SEPERATOR) 
				cRow.append(__OUTPUT_LOW_COVERAGE_VALUE) #MinMaxFirstVarf
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append(__OUTPUT_LOW_COVERAGE_VALUE) #MinMaxSecondVarf
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append(__OUTPUT_LOW_COVERAGE_VALUE) #FirstLastChangeVarf
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				cRow.append(__OUTPUT_LOW_COVERAGE_VALUE) #MinMaxChangeVarf
				cRow.append(__OUTPUT_ROW_SEPERATOR)
			else:
				if(varPercMin == 0):
					varPercMin_PlusOne = varPercMin + 1
					varPercMax_PlusOne = varPercMax + 1
					varPercMinMaxChange = ((varPercMax_PlusOne - varPercMin_PlusOne) / varPercMin_PlusOne) * 100.0 #not(varPercMax == 0) # Since 0 to x% is always >100%
				else:
					varPercMinMaxChange = ((varPercMax - varPercMin) / varPercMin) * 100.0
				IsPotentialMajority = varPercMax >= Minimal_Value_Row_Var_Perc
				IsChangeExceeding = abs(varPercMinMaxChange) >= Minimal_Value_Percentage_Change
				IsPercDifferenceExceeding = abs(varPercMax - varPercMin) >= Minimal_Value_Percentage_Difference
				IsPercDifferenceLowRegion = (abs(varPercMax - varPercMin) < Minimal_Value_Percentage_Difference) and (abs(varPercMax - varPercMin) >= Threshold_Value_Percentage_Difference)
				IsThresholdValueMaxReached = varPercMax >= Threshold_Value_Max_Reached
				
				
				# PotentialMaj
				if(IsPotentialMajority):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
					
				# PercChngOK
				if(IsChangeExceeding):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				# PercDiffOk
				if(IsPercDifferenceExceeding):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				# PercDiffLowRegion
				if(IsPercDifferenceLowRegion):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				# MaxAtLeast
				if(IsThresholdValueMaxReached):
					cRow.append("true")
				else:
					cRow.append("false")
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				

				# MinMaxFirstVarf, MinMaxSecondVarf
				bShowGreaterThan = False # MinMaxFirstVarf, MinMaxSecondVarf
				varPercChangeMinMax = 0 # MinMaxFirstVarf, MinMaxSecondVarf
				if(varPercMinSetLast):
					cRow.append(str(locale.format('%f', varPercMax))) # MinMaxFirstVarf
					cRow.append(__OUTPUT_ROW_SEPERATOR)
					cRow.append(str(locale.format('%f', varPercMin))) # MinMaxSecondVarf
					cRow.append(__OUTPUT_ROW_SEPERATOR)
					
					if(varPercMax == 0 and varPercMax == 0):
						varPercChangeMinMax = 0
					elif(varPercMax == 0):
						varPercChangeMinMax = (((varPercMin + 1) - 1) / 1) * 100.0 # +1 So we don't divide by 0
						bShowGreaterThan = True
					elif(varPercMin == 0):
						varPercChangeMinMax = -100
						bShowGreaterThan = False
					else:
						varPercChangeMinMax = ((varPercMin - varPercMax) / varPercMax) * 100.0
				else:
					cRow.append(str(locale.format('%f', varPercMin))) # MinMaxFirstVarf
					cRow.append(__OUTPUT_ROW_SEPERATOR)
					cRow.append(str(locale.format('%f', varPercMax))) # MinMaxSecondVarf
					cRow.append(__OUTPUT_ROW_SEPERATOR)
					
					if(varPercMax == 0 and varPercMax == 0):
						varPercChangeMinMax = 0
					elif(varPercMin == 0):
						varPercChangeMinMax = (((varPercMax + 1) - 1) / 1) * 100.0 # +1 So we don't divide by 0
						bShowGreaterThan = True
					elif(varPercMax == 0):
						varPercChangeMinMax = -100
						bShowGreaterThan = False
					else:
						varPercChangeMinMax = ((varPercMax - varPercMin) / varPercMin) * 100.0
					
		
				# FirstLastChangeVarf	
				if(varPercChangeFirstLast_First == varPercChangeFirstLast_First_Init or varPercChangeFirstLast_Last == varPercChangeFirstLast_Last_Init):
					cRow.append(__OUTPUT_LOW_COVERAGE_VALUE)
				elif(varPercChangeFirstLast_First == 0 and varPercChangeFirstLast_Last == 0):
					cRow.append("0")
				elif(varPercChangeFirstLast_First == 0):				
					varPercChangeFirstLast = (((varPercChangeFirstLast_Last + 1) - 1) / 1) * 100.0 # +1 So we don't divide by 0
					cRow.append(">" + str(locale.format('%f', varPercChangeFirstLast)))
				else:
					if(varPercChangeFirstLast_Last == 0):
						cRow.append("-100")
					else:
						varPercChangeFirstLast = ((varPercChangeFirstLast_Last - varPercChangeFirstLast_First) / varPercChangeFirstLast_First) * 100.0
						cRow.append(str(locale.format('%f', varPercChangeFirstLast)))
				cRow.append(__OUTPUT_ROW_SEPERATOR)
				
				
				# MinMaxChangeVarf
				if(bShowGreaterThan):
					cRow.append(">" + str(locale.format('%f', varPercChangeMinMax)))
				else:
					cRow.append(str(locale.format('%f', varPercChangeMinMax)))
				#cRow.append(__OUTPUT_ROW_SEPERATOR)
				
			
			if(not NonPrimersOnly):
				retVal.append(cRow)
			elif (not isPrimerPosition):
				retVal.append(cRow)
				
	return retVal

# Save CSV data to file	
def saveCSV(CSVTable, currentSegment, FolderName, FileNameAppendix = ""):
	file = open(outputFolder + "/" + FolderName + "/" + logfile_Name_Prefix + "_" + currentSegment + "_count" + FileNameAppendix + ".csv", "w")
	file.write("sep=" + __OUTPUT_ROW_SEPERATOR + "\n")
	
	for line in CSVTable:
		sLine = ""
		for item in line:
			sLine += str(item)
		file.write(sLine + "\n")
	
	file.close()

	

# Read the time from provided files and return a list of these files
# Provide file as None is allowed, this value will be skipped.
def GetTimePointsFromFiles(file, similarFiles):
	#Combine files into 1 list
	files = []
	if(not(file is None)):
		files.append(file)
	files.extend(similarFiles)
	
	retVal = []
	for f in files:
		cTime = GetUniqueSortingValueFromFileName(f)
		retVal.append(cTime)
	

	
	# Sort, so timelines are correct (assume time is 1st number in file name)
	retVal = sorted(retVal, key=lambda timepoint: [timepoint])
	
	# We need to change the locations values to numbers
	if(not(IsLocationInFilesUnique)):
		count = 1
		for i in range(len(retVal)):
			retVal[i] = count
			count = count + 1
	
	return retVal
	
# Process found files
def ProcessFiles(Files, PrimerLocations):
	if(len(Files) == 0):
		Show ("No files to process. Quit program")
		sys.exit(1)
		
	ProcessedFiles = []
	Files = sorted(Files, key=lambda filename: [GetUniqueSortingValueFromFileName(filename)])
	numberOfSegmentsAnalysed = 0
	
	# Get the highest timepoint over all files, so gather_positions_of_interest is nicely aligned
	highestTimePointSeen = -1
	lowestTimePointSeen = -1
	timePoints = GetTimePointsFromFiles(None, Files)

	highestTimePointSeen = max(timePoints)
	lowestTimePointSeen = min(timePoints)
	

	
	for file in Files:
		if (file in ProcessedFiles):
			continue
		# don't use this file again
		ProcessedFiles.append(file)
		
		# Find similar files
		similarFiles = FindSimilarFiles(Files, file)
		ProcessedFiles.extend(similarFiles)
		
		# Copy the correct reading frame file or place a warning file
		
		currentSegment = GetFileSegment(file)
		numberOfSegmentsAnalysed = numberOfSegmentsAnalysed + 1
		Show("Analysing %s: %s" % (DataTypeFound_DisplayName, currentSegment))
		
		timeTable = readDataFromFile(file, similarFiles, __READ_DATA_TYPE_TIME)
		readingDepth = readDataFromFile(file, similarFiles, __READ_DATA_TYPE_DEPTH)
		TimePoints = GetTimePointsFromFiles(file, similarFiles)
		CSVTable = transformTimeTableToCSV(timeTable, 1 + len(similarFiles), readingDepth, len(similarFiles) + 1, lowestTimePointSeen, highestTimePointSeen, currentSegment, PrimerLocations, TimePoints)
		saveCSV(CSVTable, currentSegment, __FILE_OUTPUT_FOLDER_WITHPRIMER)
		
		CSVTableNonPrimersOnly = transformTimeTableToCSV(timeTable, 1 + len(similarFiles), readingDepth, len(similarFiles) + 1, lowestTimePointSeen, highestTimePointSeen, currentSegment, PrimerLocations, TimePoints, True)
		saveCSV(CSVTableNonPrimersOnly, currentSegment, __FILE_OUTPUT_FOLDER_NONPRIMER, __FILE_OUTPUT_APPENDIX_NONPRIMERSONLY)
		
		#DetermineCorrectReadingFrame(file, similarFiles)
	
	Show("Analysed %s %s(s)" % (numberOfSegmentsAnalysed, DataTypeFound_DisplayName))
		
# Write log file
def WriteLogFile(Files):	
	file = open(outputFolder + "/" + __FILE_OUTPUT_FOLDER_LOG + "/" + logfile_Name_Prefix + ".log", "w")
	# file.write("Log file for segment over time tool\n")  -- BK adapted 06MAY17
	file.write("Log file for variant_detector.py tool\n")
	file.write("\n")
	file.write("SETTINGS\n")
	file.write("showMessages: " + ("TRUE" if showMessages else "FALSE")+ "\n")
	file.write("deepDebug: " + ("TRUE" if deepDebug else "FALSE") + "\n")
	file.write("DataTypeFound: " + DataTypeFound + "\n")
	file.write("tabularDivide: " + tabularDivide + "\n") 
	file.write("stopCodon: " + stopCodon + "\n")
	file.write("Minimal_Value_RowCoverage: " + str(Minimal_Value_RowCoverage) + "\n") 
	file.write("__VarCovOK_Minimal_Value_Count: " + str(__VarCovOK_Minimal_Value_Count) + "\n") 
	file.write("__VarCovOK_Minimal_Value_Percentage: " + str(__VarCovOK_Minimal_Value_Percentage) + "\n") 
	file.write("Minimal_Value_Percentage_Change: " + str(Minimal_Value_Percentage_Change) + "\n") 
	file.write("Minimal_Value_Percentage_Difference: " + str(Minimal_Value_Percentage_Difference) + "\n") 
	file.write("Minimal_Value_Row_Var_Perc: " + str(Minimal_Value_Row_Var_Perc) + "\n") 
	
	file.write("Threshold_Value_Percentage_Difference: " + str(Threshold_Value_Percentage_Difference) + "\n") 
	file.write("Threshold_Value_Max_Reached: " + str(Threshold_Value_Max_Reached) + "\n") 
	
	file.write("primerFile: " + primerFile + "\n")
	file.write("inputFolder: " + inputFolder + "\n")
	file.write("outputFolder: " + outputFolder + "\n")
	file.write("logfile_Name_Prefix: " + logfile_Name_Prefix + "\n")
	file.write("__FILE_CONTENT_COL_ROW: " + str(__FILE_CONTENT_COL_ROW) + "\n")
	file.write("__FILE_CONTENT_COL_DEPTH: " + str(__FILE_CONTENT_COL_DEPTH)+ "\n")
	file.write("__OUTPUT_ROW_SEPERATOR: " + __OUTPUT_ROW_SEPERATOR+ "\n")
	
	file.write("__DATA_LETTERS_SELECTED: ")
	for dt in __DATA_LETTERS_SELECTED: 
		file.write(dt + " ")
	file.write("\n")

	file.write("\n")
	file.write("FILES\n")
	
	
	Files = sorted(Files, key=lambda filename: [GetUniqueSortingValueFromFileName(filename)])
	
	count = 1
	for f in Files:
		file.write("f" + str(count) + ": " + f + "\n")
		count = count + 1

	file.close()

def GetPrimerLocations():
	primerLocations = {}
	
	for line in open(primerFile):
		data = line.split(";")
		position = int(data[1].rstrip())
		
		if(not (data[0] in primerLocations)):
			primerLocations[data[0]] = [position]
		else:
			primerLocations[data[0]].append(position)
	
	return primerLocations

# Initialise arguments
def Initialise():
	# Write numbers with "," instead of "."
	# locale.setlocale(locale.LC_ALL, 'nl_NL')
	
	args = len(sys.argv)
	if(not (args >= 5)):
		QuitAndShowError([("Not all arguments are provided")])
		
	global primerFile
	global inputFolder
	global outputFolder
	global logfile_Name_Prefix
	
	global Minimal_Value_RowCoverage
	global __VarCovOK_Minimal_Value_Count
	global __VarCovOK_Minimal_Value_Percentage
	global Minimal_Value_Percentage_Change
	global Minimal_Value_Row_Var_Perc
	global Minimal_Value_Percentage_Difference
	global Threshold_Value_Percentage_Difference
	global Threshold_Value_Max_Reached
	global __RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time
	global __RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location

	inputFolder = sys.argv[1]
	outputFolder = sys.argv[2]
	primerFile = sys.argv[3]
	logfile_Name_Prefix = sys.argv[4] 

	if(args > 5):
		Show("Using configuration send by parameters.")
		Minimal_Value_RowCoverage = int(sys.argv[5])
		__VarCovOK_Minimal_Value_Count = int(sys.argv[6])
		__VarCovOK_Minimal_Value_Percentage = int(sys.argv[7])
		Minimal_Value_Percentage_Change = int(sys.argv[8])
		Minimal_Value_Row_Var_Perc = int(sys.argv[9])
		Minimal_Value_Percentage_Difference = int(sys.argv[10])
		Threshold_Value_Percentage_Difference = int(sys.argv[11])
		Threshold_Value_Max_Reached = int(sys.argv[12])
		__RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time = int(sys.argv[13])
		__RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location = int(sys.argv[14])
	else:
		Show("No configuration parameters provided. Continue with default values specified in this python script.")

	
def InitialiseDataTypeLetters():
	global DataTypeFound
	global __DATA_LETTERS_SELECTED
	global __FIRST_COUNT_COL_HEADER_NAME_SELECTED
	
	# AminoAcids or Nucleotides?
	if (DataTypeFound == __DATATYPE_AMINOACIDS):
		__DATA_LETTERS_SELECTED = __DATA_LETTERS_AA
		__FIRST_COUNT_COL_HEADER_NAME_SELECTED = __FIRST_COUNT_COL_HEADER_NAME_AA
	elif (DataTypeFound == __DATATYPE_NUCLEOTIDES):
		__DATA_LETTERS_SELECTED = __DATA_LETTERS_NT
		__FIRST_COUNT_COL_HEADER_NAME_SELECTED = __FIRST_COUNT_COL_HEADER_NAME_NT
	else:
		QuitAndShowError([("Detected datatype is not supported."),
						("Detected: " + DataTypeFound),
						("Supported: [" + __DATATYPE_AMINOACIDS + "/" + __DATATYPE_NUCLEOTIDES + "]")])
						
	
def MakeOutputFolders():
	Show("Make folders:")
	Show("\t-" + outputFolder + "/" + __FILE_OUTPUT_FOLDER_LOG)
	if not os.path.exists(outputFolder + "/" + __FILE_OUTPUT_FOLDER_LOG):
		os.makedirs(outputFolder + "/" + __FILE_OUTPUT_FOLDER_LOG)
		
	Show("\t-" + outputFolder + "/" + __FILE_OUTPUT_FOLDER_NONPRIMER)
	if not os.path.exists(outputFolder + "/" + __FILE_OUTPUT_FOLDER_NONPRIMER):
		os.makedirs(outputFolder + "/" + __FILE_OUTPUT_FOLDER_NONPRIMER)
		
	Show("\t-" + outputFolder + "/" + __FILE_OUTPUT_FOLDER_WITHPRIMER)
	if not os.path.exists(outputFolder + "/" + __FILE_OUTPUT_FOLDER_WITHPRIMER):
		os.makedirs(outputFolder + "/" + __FILE_OUTPUT_FOLDER_WITHPRIMER)
		
	Show("\t-" + outputFolder + "/" + __FILE_OUTPUT_FOLDER_CONSENSUS)
	if not os.path.exists(outputFolder + "/" + __FILE_OUTPUT_FOLDER_CONSENSUS):
		os.makedirs(outputFolder + "/" + __FILE_OUTPUT_FOLDER_CONSENSUS)

# Check if files are all the same times and equal locations or visa versa
def AutoDetectIfTimeOrLocationIsUnique(Files):
	global IsTimeInFilesUnique
	global IsLocationInFilesUnique
	global FileNames_UniqueValues
	
	uniqueValues_Time = []
	uniqueValues_Location = []
	
	TimeIsDifferent = False
	FirstTimeFound = -1
	for file in Files:
		cTime = GetUniqueSortingValueFromFileName(file,__FILE_NAME_TIME_COL, __FILE_NAME_TIME_PRE_SIZE)
		uniqueValues_Time.append(cTime)
		if (FirstTimeFound == -1):
			FirstTimeFound = cTime
		else:
			if (not(cTime == FirstTimeFound)):
				TimeIsDifferent = True
	
	# Check locations
	LocationIsDifferent = False
	FirstLocationFound = ""
	for file in Files:
		cLocation = GetUniqueSortingValueFromFileName(file ,__FILE_NAME_LOCATION_COL, __FILE_NAME_LOCATION_PRE_SIZE, False)
		uniqueValues_Location.append(cLocation)
		if (FirstLocationFound == ""):
			FirstLocationFound = cLocation
		else:
			if (not(cLocation == FirstLocationFound)):
				LocationIsDifferent = True
			
	
	# If we did find different time markers, do not accept different locations	
	if (TimeIsDifferent and LocationIsDifferent):
		err = [("(CheckIfFilesAreLegal) Provided data have different time and location markers. This is not supported in this current version."),
							("Files provided:")]
		for file in Files:
			err.append(" - " + file)
		QuitAndShowError(err)

	
	if(not(TimeIsDifferent)):
		IsTimeInFilesUnique = True
		IsLocationInFilesUnique = False
		FileNames_UniqueValues = uniqueValues_Location
	else:
		IsTimeInFilesUnique = False
		IsLocationInFilesUnique = True
		FileNames_UniqueValues = uniqueValues_Time
	
	FileNames_UniqueValues = set(FileNames_UniqueValues)
	FileNames_UniqueValues = sorted(FileNames_UniqueValues, key=lambda uniqueValue: [uniqueValue])
	
		
def getConsensus(line):
	startConsensus = find_nth(line,tabularDivide, __FILE_CONTENT_COL_CONSENSUS)
	endConsensus = find_nth(line,tabularDivide, __FILE_CONTENT_COL_CONSENSUS + 1)
	return line[startConsensus + 1: endConsensus]
			
def GetConsensusSequence(File, ReadingDepthThreshold = 0):
	consensus = ""
	isHeader = True
	hasConsensusAboveThreshold = False
	for line in open(inputFolder + "/" + File):
		# Skip Header line
		if(isHeader):
			isHeader = False
			continue
		depth = GetDepthFromFileLine(line)
		
		if (depth >= ReadingDepthThreshold):
			consensus = consensus + getConsensus(line)
			hasConsensusAboveThreshold = True
		else:
			consensus = consensus + "x"
	return consensus, hasConsensusAboveThreshold
	
def BuildConsensusToFasta(Files):
	fileFastaContent = ""
	
	outputFile = outputFolder + "/" + __FILE_OUTPUT_FOLDER_CONSENSUS + "/" + __FILE_OUTPUT_FILENAME_CONSENSUS
	Show ("Build consensus overview file: " + outputFile)
	
	for file in Files:
		consensus, hasConsensusAboveThreshold = GetConsensusSequence(file, Minimal_Value_RowCoverage)
		if (hasConsensusAboveThreshold):
			fileFastaContent = fileFastaContent + ">" + file + "\n"
			fileFastaContent = fileFastaContent + consensus + "\n"
	
	
	
	file = open(outputFile, "w")
	file.write(fileFastaContent)
	file.close()
	
# Main method
def Main():

	ShowSplash()
	Initialise()
	MakeOutputFolders()
	
	Files = GetFileDirectPaths(inputFolder)
	Files = RemoveNonReadingFramedFiles(Files)
	AutoDetectIfTimeOrLocationIsUnique(Files)

	BuildConsensusToFasta(Files)
	PrimerLocations = GetPrimerLocations()


	InitialiseDataTypeLetters()
	WriteLogFile(Files)
	ProcessFiles(Files, PrimerLocations)
	
	Show("Ready!") 

Main()