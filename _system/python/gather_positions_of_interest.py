#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
#  Gather positions of interests		                    #
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
		
# Constants
showMessages = True
deepDebug = True
DataTypeFound = ""
ChangeNonValue = "-"
__DATATYPE_AMINOACIDS = "AA"
__DATATYPE_NUCLEOTIDES = "NA"
tabularDivide =  "\t"
		
inputFolderNT = ""
inputFolderAA = ""
outputFolder = ""
IdentifiableExperimentName = ""
		
__FILE_NAME_TIME_DIVIDE = "_" # used to get time from filename
		
__FILE_OUTPUT_APPENDIX_NONPRIMERSONLY = ".NonPrimersOnly"
__FILE_OUTPUT_APPENDIX_NT = ".NT"
__FILE_OUTPUT_APPENDIX_AA = ".AA"
__FILE_CONTENT_COL_SEGMENT = 1
# __FILE_CONTENT_COL_COVERAGE_OK = 6
# __FILE_CONTENT_COL_VAR_COVERAGE_OK = 7
# __FILE_CONTENT_COL_PRIMER_POSITION = 8
__FILE_CONTENT_COL_POSITION = 2
# __FILE_CONTENT_HEADER_COVERAGE_ROW_OK = "RowCoverageOK"
__FILE_CONTENT_HEADER_VAR_COVERAGE_OK = "RowMultiVarCovOK"
__FILE_CONTENT_HEADER_POTENTIALMAJORITY_OK = "PotentialMaj"
__FILE_CONTENT_HEADER_PERCCHNG_OK = "PercChngOK"
__FILE_CONTENT_HEADER_PRIMER_POSITION = "IsPrimer"
# __FILE_CONTENT_HEADER_CHANGE_STARTS_WITH = "VarPercChange"
__FILE_SEGMENTED_OVERVIEW_FILE_EXTENTION = "csv"
__FILE_SEGMENTED_OVERVIEW_FILE_NAME_IDENTIFIABLE_END_CHAR = "-"
		
		
# Show some startup information
def ShowSplash():
	Show ("")
	Show ("--------------------------------")
	Show ("Tool: Gather positions of interest")
	Show ("Read all segmented overviews CSV files and gather only interesting rows.")
	Show ("version: 1.0")
	Show ("Contact: m.pater@amc.nl | mpater@mirdesign.nl - Maarten Pater")
	Show ("--------------------------------")
	Show ("Usage: python gather_positions_of_interests.py folder_with_segmented_overview_files_NT folder_with_segmented_overview_files_AA output_folder experiment_name")
	Show ("Please provide the primer_incl folders from variant_detector.")
	Show ("--------------------------------")
		
			
# Detect header locations
def GetColumnPositionFromHeader(line, SearchForHeader):
	pos = line.find(SearchForHeader)
			
	# Check if we found the headers
	if (pos == -1):
		QuitAndShowError([	("(GetColumnPositionFromHeader) " + SearchForHeader + " value could not be detected in the headerline."),
							("Provided: " + line),
							("Searched for: " + SearchForHeader)])
							
	return line[:pos].count(__INPUT_ROW_SEPERATOR)
	
									
# Transform timetable from file to count output	
def ExtractPositionsOfInterests(file, inputFolder, NonPrimersOnly = False):
	retVal = []
	potentialMajority = []
	
	skippedFirstLine = False
	headerLine = True
				
	# colRowCoverageOK = -1
	colRowMultiVarCoverageOK = -1
	colPercChngOK = -1
	colPrimer = -1
	
	isPositionInteresting = False
	lastPosition = -1
	
	for line in open(inputFolder + "/" + file):
		# Skip seperation notation (Propably: "Sep=;" used for Excel)
		if(not skippedFirstLine):
			skippedFirstLine = True
			continue
					
		if(headerLine):
			# Get header position from header
			# colRowCoverageOK = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_COVERAGE_ROW_OK)
			colRowMultiVarCoverageOK = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_VAR_COVERAGE_OK)
			colPotMaj = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_POTENTIALMAJORITY_OK)
			colPercChngOK = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_PERCCHNG_OK)
			
			if(NonPrimersOnly):
				colPrimer = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_PRIMER_POSITION)
					
			headerLine = False
			continue
				
				
		currentPosition = int(GetDataFromFileContent(line, __FILE_CONTENT_COL_POSITION))
		if (lastPosition != currentPosition):
			# If we found potential majorities on the interesting position, add them to the results
			if(isPositionInteresting):
				if(len(potentialMajority) > 0):
					retVal.extend(potentialMajority)
			isPositionInteresting = False
			potentialMajority = []
		
		# Get isCoverageOk value
		# isRowCoverageOKValue = GetDataFromFileContent(line, colRowCoverageOK).lower()
		isRowMultiVarCoverageOKValue = GetDataFromFileContent(line, colRowMultiVarCoverageOK).lower()
		isPercChngOK = GetDataFromFileContent(line, colPercChngOK).lower()
		isPrimerPosition = GetDataFromFileContent(line, colPrimer).lower()
		isPotentialMajority = GetDataFromFileContent(line, colPotMaj).lower()
		# if(isPercChngOK == "true" and isRowCoverageOKValue == "true" and isRowMultiVarCoverageOKValue == "true"):
		if(isPercChngOK == "true" and isRowMultiVarCoverageOKValue == "true"):
			if((not NonPrimersOnly) or (NonPrimersOnly and (not isPrimerPosition == "true"))):
				if (line[-1:] == "\n"):
					retVal.append(line[0:len(line) - 1])
					isPositionInteresting = True
				else:
					retVal.append(line)
					isPositionInteresting = True
		else:
			# Check if this Data might be interesting if we find a interesting datapoint on the same position
			# So remember, and append to the output if we found an interesting position
			if(isPotentialMajority == "true"):
				if (line[-1:] == "\n"):
					potentialMajority.append(line[0:len(line) - 1])
				else:
					potentialMajority.append(line)
				
			
			
		lastPosition = int(GetDataFromFileContent(line, __FILE_CONTENT_COL_POSITION))
			
	if(isPositionInteresting):
		if(len(potentialMajority) > 0):
			retVal.extend(potentialMajority)
		
	return retVal, lastPosition
		
		
# Save CSV data to file	
def saveCSV(CSVTable, FileNameAppendix = ""):
	file = open(outputFolder + "/" + IdentifiableExperimentName + FileNameAppendix + ".csv", "w")
	file.write("sep=" + __OUTPUT_ROW_SEPERATOR + "\n")
	for line in CSVTable:
		sLine = ""
		for item in line:
			sLine += str(item)
		file.write(sLine + "\n")
			
	file.close()
		
def CountNumberOfInterestingPositions(currentPOIs):
	lastCol = -1
	numberOfPOIPositions = 0
	for line in currentPOIs:
		currentCol = int(GetDataFromFileContent(line, __FILE_CONTENT_COL_POSITION))
		if (not currentCol == lastCol):
			lastCol = currentCol
			numberOfPOIPositions = numberOfPOIPositions + 1
					
	return numberOfPOIPositions
		
def getHeaderLine(file, inputFolder):
	skipFirstLine = True
	retVal = []
	for line in open(inputFolder + "/" + file):
		if (skipFirstLine):
			skipFirstLine = False
			continue
		
		if (line[-1:] == "\n"):
			retVal.append(line[0:len(line) - 1])
		else:
			retVal.append(line)
		return retVal
	return ["Header information not found"]
				
# Process found files
def ProcessFiles(Files, inputFolder, NonPrimersOnly = False):
	if(len(Files) == 0):
		Show ("No files to process. Quit program")
		sys.exit(1)
			
	POIData = []
	overviewData = []
			
			
	for file in Files:
		# Copy the correct reading frame file or place a warning file
				
		Show("Analysing file: %s" % file)
		currentPOIs, lastPosition = ExtractPositionsOfInterests(file, inputFolder, NonPrimersOnly)
				
		if (len(currentPOIs) > 0):
			POIData.extend(currentPOIs)
			currentSegment = GetDataFromFileContent(currentPOIs[0], __FILE_CONTENT_COL_SEGMENT)
					
			# Get number of unique positions
			numberOfInterestingPositions = CountNumberOfInterestingPositions(currentPOIs)
			percentagePOIs = (100.0 / lastPosition) * numberOfInterestingPositions
					
			overviewDataLine = {}
			overviewDataLine["segment"] = currentSegment
			overviewDataLine["numberOfInterestingPositions"] = (str(numberOfInterestingPositions))
			overviewDataLine["percentagePOIs"] = (str(percentagePOIs))
					
			overviewData.append(overviewDataLine)	
		
		
			
	return POIData, overviewData
			
# Write log file
def WriteLogFile(FilesNT, FilesAA):	
		
	file = open(outputFolder + "/" + IdentifiableExperimentName + ".log", "w")
	file.write("Log file for gather positions of interests tool\n")
	file.write("\n")
	file.write("SETTINGS\n")
	file.write("showMessages: " + ("TRUE" if showMessages else "FALSE")+ "\n")
	file.write("deepDebug: " + ("TRUE" if deepDebug else "FALSE") + "\n")
	file.write("DataTypeFound: " + DataTypeFound + "\n")
	file.write("tabularDivide: " + tabularDivide + "\n") 
			
	file.write("inputFolder NT: " + inputFolderNT + "\n")
	file.write("inputFolder AA: " + inputFolderAA + "\n")
	file.write("outputFolder: " + outputFolder + "\n")
	file.write("IdentifiableExperimentName: " + IdentifiableExperimentName + "\n")
	file.write("__FILE_CONTENT_COL_POSITION: " + str(__FILE_CONTENT_COL_POSITION) + "\n")
	file.write("__FILE_SEGMENTED_OVERVIEW_FILE_EXTENTION: " + str(__FILE_SEGMENTED_OVERVIEW_FILE_EXTENTION)+ "\n")
	file.write("__OUTPUT_ROW_SEPERATOR: " + str(__OUTPUT_ROW_SEPERATOR)+ "\n")
			
			
	file.write("\n")
	file.write("FILES\n")
	for f in FilesNT:
		file.write(inputFolderNT + "/" + f + "\n")
	for f in FilesAA:
		file.write(inputFolderAA + "/" + f + "\n")
		
	file.close()
		
# Initialise arguments
def Initialise():
	# Write numbers with "," instead of "."
	# locale.setlocale(locale.LC_ALL, 'nl_NL')
			
	args = len(sys.argv)
	if(not (args >= 5)):
		Show("ERROR Not all arguments are provided")
		sys.exit(1)
				
	global inputFolderAA
	global inputFolderNT
	global outputFolder
	global IdentifiableExperimentName
			
	inputFolderNT = sys.argv[1]
	inputFolderAA = sys.argv[2]
	outputFolder = sys.argv[3]
	IdentifiableExperimentName = sys.argv[4]
	
		
# Removing any file which isn't AA/nt typed
def RemoveNonSegmentedOverviewFiles(Files):
			
	global DataTypeFound
	remainingFiles = []
				
	for file in Files:
		if((file[-3:] == __FILE_SEGMENTED_OVERVIEW_FILE_EXTENTION)): # Skip log files, probably created by reading_frame_selector.py
			remainingFiles.append(file)
			
	Show("Found %s typed files: %s" % (__FILE_SEGMENTED_OVERVIEW_FILE_EXTENTION, str(len(remainingFiles))))
	return remainingFiles
		
def WriteOverviewData(overviewDataNT, overviewDataAA, FileNameAppendix = ""):
	file = open(outputFolder + "/" + IdentifiableExperimentName + "_experiment_overview" + FileNameAppendix + ".csv", "w")
	file.write("sep=" + __OUTPUT_ROW_SEPERATOR + "\n")
	file.write("segment")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("NT #pos")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("AA #pos")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("NT %pos")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("AA %pos")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("\n")
			
	foundSegments = []
	for segmentDataNT in overviewDataNT:
		foundSegments.append(segmentDataNT["segment"])
	for segmentDataAA in overviewDataAA:
		if(not segmentDataAA["segment"] in foundSegments):
			foundSegments.append(segmentDataAA["segment"])

	for segment in foundSegments:
		#segment = segmentDataNT["segment"]
		segmentDataNT = None
		segmentDataAA = None
		for segmentDataNTItem in overviewDataNT:
			if (segmentDataNTItem["segment"] == segment):
				segmentDataNT = segmentDataNTItem
				break
				
		for segmentDataAAItem in overviewDataAA:
			if (segmentDataAAItem["segment"] == segment):
				segmentDataAA = segmentDataAAItem
				break
			
		file.write(segment)
		file.write(__OUTPUT_ROW_SEPERATOR)
		if(segmentDataNT == None):
			file.write("0")
		else:
			file.write(segmentDataNT["numberOfInterestingPositions"])
		file.write(__OUTPUT_ROW_SEPERATOR)
		if(segmentDataAA == None):
			file.write("0")
		else:
			file.write(segmentDataAA["numberOfInterestingPositions"])
		file.write(__OUTPUT_ROW_SEPERATOR)
		if(segmentDataNT == None):
			file.write("0")
		else:
			file.write(segmentDataNT["percentagePOIs"])
		file.write(__OUTPUT_ROW_SEPERATOR)
				
		if(segmentDataAA == None):
			file.write("0")
		else:
			file.write(segmentDataAA["percentagePOIs"])
		file.write(__OUTPUT_ROW_SEPERATOR)
		file.write("\n")
				
				
def GetIdentifiableNameFromFileName(fileName):
	endOfIdentifiableRegion = fileName.find(__FILE_SEGMENTED_OVERVIEW_FILE_NAME_IDENTIFIABLE_END_CHAR)
	return fileName[0:endOfIdentifiableRegion]
		
def CheckIfAllFilesAreFromTheSameExperiment(FilesNT, FilesAA):
	FilesToCheck = []
	FilesToCheck.extend(FilesNT)
	FilesToCheck.extend(FilesAA)
			
	IdentifiableNameBase = GetIdentifiableNameFromFileName(FilesToCheck[0])
			
	for fileToCheck in FilesToCheck:
		IdentifiableNameCurrent = GetIdentifiableNameFromFileName(fileToCheck)
		if(not IdentifiableNameBase == IdentifiableNameCurrent):
			Show("ERROR (CheckIfAllFilesAreTheSameExperiment) At least 2 experimental output files are provided.")
			Show("ERROR IdentifiableNameBase: " + str(IdentifiableNameBase))
			Show("ERROR IdentifiableNameCurrent: " + str(IdentifiableNameCurrent))
			exit(1)
		
# Main method
def Main():
		
	ShowSplash()
	Initialise()
			
	# NT
	Show("")
	Show("Searching for NT files...")
	FilesNT = GetFileDirectPaths(inputFolderNT)
	FilesNT = RemoveNonSegmentedOverviewFiles(FilesNT)
			
	# AA
	Show("")
	Show("Searching for AA files...")
	FilesAA = GetFileDirectPaths(inputFolderAA)
	FilesAA = RemoveNonSegmentedOverviewFiles(FilesAA)

	#######################################
	# Gather POI positions - incl primers #
	#######################################
	Show("")
	Show ("Gather positions of interests, including primer positions...")
	POIDataNT, overviewDataNT = ProcessFiles(FilesNT, inputFolderNT)
	POIDataAA, overviewDataAA = ProcessFiles(FilesAA, inputFolderAA)
	WriteOverviewData(overviewDataNT, overviewDataAA)
		
	Show("Write output file...")
	WriteLogFile(FilesNT, FilesAA)
			
	# Write POI positions to 2 files
	# NT
	allPOIDataNT = [getHeaderLine(FilesNT[0], inputFolderNT)]
	allPOIDataNT.extend(POIDataNT)
	saveCSV(allPOIDataNT, __FILE_OUTPUT_APPENDIX_NT)
			
	# AA
	allPOIDataAA = [getHeaderLine(FilesAA[0], inputFolderAA)]
	allPOIDataAA.extend(POIDataAA)
	saveCSV(allPOIDataAA, __FILE_OUTPUT_APPENDIX_AA)
			
			
			
	########################################################
	# Gather POI positions - with Primer locations removed #
	########################################################
	Show("")
	Show ("Gather positions of interests, excluding primer positions...")
	POIDataNT, overviewDataNT = ProcessFiles(FilesNT, inputFolderNT, True)
	POIDataAA, overviewDataAA = ProcessFiles(FilesAA, inputFolderAA, True)
	WriteOverviewData(overviewDataNT, overviewDataAA, __FILE_OUTPUT_APPENDIX_NONPRIMERSONLY)
			
	Show("Write output file...")
			
	# Write POI positions to 2 files
	# NT
	allPOIDataNT = [getHeaderLine(FilesNT[0], inputFolderNT)]
	allPOIDataNT.extend(POIDataNT)
	saveCSV(allPOIDataNT, __FILE_OUTPUT_APPENDIX_NT + __FILE_OUTPUT_APPENDIX_NONPRIMERSONLY)
			
	# AA
	allPOIDataAA = [getHeaderLine(FilesAA[0], inputFolderAA)]
	allPOIDataAA.extend(POIDataAA)
	saveCSV(allPOIDataAA, __FILE_OUTPUT_APPENDIX_AA + __FILE_OUTPUT_APPENDIX_NONPRIMERSONLY)
				
				
	Show("Ready!") 
		
Main()