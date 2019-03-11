#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# Protein numbering convertor                               #
#############################################################
from __future__ import print_function
import sys 
import os
import glob
from shutil import copyfile

# Load common shared variables and methods
execfile(os.path.join(os.path.dirname(__file__),"_global_constants.py"))

inputFolder = ""
inputFolderNT = ""
outputFolder = ""
outputFolderNT = ""
outputFolderNT = ""
outputWarningLogFile = "warnings.log"
_warning_file_handler = ""

COLUMNS_POSITION = 1 # Also 1, as old Position, since this has been added by the same script
COLUMNS_CONSENSUS = 2
COLUMNS_READINGDEPTH = 3
COLUMNS_ADD_TO_POSITION = COLUMNS_READINGDEPTH + 1 # COLUMNS_ADD_TO_POSITION: Always behind columns we add ourself
HEADER_PROTEINCOUNT = "P_Cnt"
HEADER_H3 = "H3"
HEADER_N40 = "N40"

SEGMENTS_DESIRED_LENGTHS = {}
PROTEIN_GAPS = {}


# Write message to warning log
def writeToWarningLog(message, showToScreen = False):
	_warning_file_handler.write(message)
	_warning_file_handler.write("\n")
	if(showToScreen):
		Show(message)

# Write array to warning log
def  writeArrayToWarningLog(messageArray, showToScreen = False):
	for message in messageArray:
		writeToWarningLog(message, showToScreen)

# Show some startup information
def ShowSplash():
	Debug ("")
	Debug ("--------------------------------")
	Debug ("Tool: Protein numbering converter")
	Debug ("Purpose: Alter influenza readingframed files and add new information")
	Debug ("version: 1.0")
	Debug ("Contact: m.pater@amc.nl | mpater@mirdesign.nl - Maarten Pater")
	Debug ("--------------------------------")
	Debug ("Usage: python protein_numbering_converter.py folder_with_reading_frame_files_AA_Type folder_with_reading_frame_files_NT_Type output_folder_AA output_folder_NT")
	Debug ("Expected Filename strategy: [__CUSTOM__]-[segment]_[" + __DATATYPE_AMINOACIDS + "/" + __DATATYPE_NUCLEOTIDES + "].ext")
	Debug ("Overwrite existing output files: " + ("Yes" if OverWriteOutputFiles else "No"))
	Debug ("--------------------------------")

	
# Get AA or nt type information
def GetFileDNAProteinType(file):
	return file[file.rfind(FILENAMES_DATATYPE_LASTDIVIDE_CHAR) + len(FILENAMES_DATATYPE_LASTDIVIDE_CHAR): file.rfind(FILENAMES_DATATYPE_LASTDIVIDE_CHAR) + len(__DATATYPE_AMINOACIDS) + len(FILENAMES_DATATYPE_LASTDIVIDE_CHAR)]

# Get position from line
def getPositionCount(line):
	start = find_nth(line,tabularDivide, COLUMNS_POSITION)
	end = find_nth(line,tabularDivide, COLUMNS_POSITION + 1)
	return int(line[start + 1: end])
	
# Get consensus value from line
def getConsensus(line):
	startConsensus = find_nth(line,tabularDivide, COLUMNS_CONSENSUS)
	endConsensus = find_nth(line,tabularDivide, COLUMNS_CONSENSUS + 1)
	return line[startConsensus + 1: endConsensus]
			
# Get consensus from file
ConsensusPerFile = {} # Cache consensus per file
def GetConsensusSequence(File):
	global ConsensusPerFile
	if(File in ConsensusPerFile):
		return ConsensusPerFile[File]

	consensus = ""
	isHeader = True
	for line in open(inputFolder + "/" + File):
		# Skip Header line
		if(isHeader):
			isHeader = False
			continue
		consensus = consensus + getConsensus(line)
	
	ConsensusPerFile[File] = consensus
	return consensus

# Get line from a file
def GetLine(File, lineNumber):
	lineCount = 0
	for line in open(inputFolder + "/" + File):
		if (lineCount == lineNumber):
			return line
		lineCount = lineCount + 1
	
	QuitAndShowError([("File does not contain the line which was requested."),
		("Searched for line: " + str(lineNumber)),
		("In file: " + File),
		("Number of lines found: " + str(lineCount))])
		
	return -1 # Will never happen, just to keep the compiler happy

# Get percentage of positions which are above minimal reading depth threshold used in variant detector (Therefor, protein is named in this function name)
# offset will be used to start from this position (ie Used for N40)
def GetPercentageAboveProteinReadingDepthThreshold(File, offset):
	offset = offset + 1 # skip header
	numberOfPositionsAboveMinimalReadingDepth = 0

	currentLine = 0
	consensusSequencesLength = 0
	for line in open(inputFolder + "/" + File):
		if currentLine >= 1 and currentLine >= offset : # Skip header and start from offset
			consensusSequencesLength = consensusSequencesLength + 1
			cReadingDepth = getReadingDepth(GetLine(File, currentLine))
			if cReadingDepth >= rejectOnStopCodonBeforeLine:
				numberOfPositionsAboveMinimalReadingDepth = numberOfPositionsAboveMinimalReadingDepth + 1

		currentLine = currentLine + 1
		
	percentageAboveMinimalReadingDepth = 0
	if(numberOfPositionsAboveMinimalReadingDepth > 0):
		percentageAboveMinimalReadingDepth = (float(numberOfPositionsAboveMinimalReadingDepth) / float(consensusSequencesLength)) * 100

	return percentageAboveMinimalReadingDepth
	
# Get the start position based on start or stop codon
def GetProteinOrSegmentCountStartPosition(File, consensusSequences):
	searchFor = ""
	
	offset = 0
	readPerLine = -9999
	segment = GetFileSegment(File)
	
	# Determine data type
	DNAProteinType = GetFileDNAProteinType(File)
	if (DNAProteinType == __DATATYPE_AMINOACIDS):
		searchFor = PROTEIN_START_AA.upper()
		readPerLine = 1
	elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
		searchFor = PROTEIN_START_NT.upper()
		readPerLine = 3
	else:
		QuitAndShowError([("File describes datatype which is not supported."),
		("Data type described as: " + DNAProteinType),
		("For file: " + File)])
	
	if(segment == SEGMENTS_PB1N40):
		offset = (SEGMENTS_PB1_SECOND_M_POS - 1) * readPerLine
	elif(segment == SEGMENTS_PB1F2):
		offset = (SEGMENTS_PB1F2_SECOND_M_POS - 1) * readPerLine
	
	foundStartPositionWithMinimalDepth = False
	startProteinCountPos = offset - readPerLine
	consensusSequences = consensusSequences.upper()
	
	while(not(foundStartPositionWithMinimalDepth)): # Will never end, but GetProteinCountStartPosition will return an error if no other startProtein can be found
		# Get start position
		if (DNAProteinType == __DATATYPE_AMINOACIDS):
			startProteinCountPos = consensusSequences.find(searchFor, startProteinCountPos + readPerLine)
		elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
			posFound = False
			overflow = False
			while(not(posFound) and not(overflow)):
				if(consensusSequences[startProteinCountPos + readPerLine: (startProteinCountPos + readPerLine) + readPerLine] == searchFor):
					posFound = True
					startProteinCountPos = startProteinCountPos + readPerLine # Just in case the readingdepth is to low. We would like to search for the next match, not the same.
				else:
					startProteinCountPos = startProteinCountPos + readPerLine
					if(startProteinCountPos > len(consensusSequences)):
						overflow = True
						startProteinCountPos = -9999
		
		if (startProteinCountPos < 0):
			startProteinCountPos = -9999
			break
		
		if (DNAProteinType == __DATATYPE_AMINOACIDS):
			readingDepth = getReadingDepth(GetLine(File, startProteinCountPos + 1)) # +1 skip header
		elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
			readingDepth1 = getReadingDepth(GetLine(File, startProteinCountPos + 1)) # +1 skip header
			readingDepth2 = getReadingDepth(GetLine(File, startProteinCountPos + 1 + 1)) # +1 skip header
			readingDepth3 = getReadingDepth(GetLine(File, startProteinCountPos + 2 + 1)) # +1 skip header
			readingDepth = min(readingDepth1, readingDepth2, readingDepth3)
		
		
		if (readingDepth < readMinimalCodonDepth):
			# Ok, we didn't found a proper M, start searching for a proper * (stop codon)
			endProteinCountPos = 0
			if (DNAProteinType == __DATATYPE_AMINOACIDS):
				endProteinCountPos = consensusSequences.find(stopCodon, endProteinCountPos + readPerLine)
			elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
				posFound = False
				overflow = False
				while(not(posFound) and not(overflow)):
					if((consensusSequences[endProteinCountPos + readPerLine: (endProteinCountPos + readPerLine) + readPerLine] == stopCodon_nt_1) or
						(consensusSequences[endProteinCountPos + readPerLine: (endProteinCountPos + readPerLine) + readPerLine] == stopCodon_nt_2) or
						(consensusSequences[endProteinCountPos + readPerLine: (endProteinCountPos + readPerLine) + readPerLine] == stopCodon_nt_3)):
						posFound = True
						endProteinCountPos = endProteinCountPos + readPerLine # Just in case the readingdepth is to low. We would like to search for the next match, not the same.
					else:
						endProteinCountPos = endProteinCountPos + readPerLine
						if(endProteinCountPos > len(consensusSequences)):
							overflow = True
							endProteinCountPos = -1
			
			# we found a *, now count back and check if that position has a readingDepth lower than our threshold
			if (endProteinCountPos == -1):
				startProteinCountPos = -9999
				break
			
			if (DNAProteinType == __DATATYPE_AMINOACIDS):
				readingDepth = getReadingDepth(GetLine(File, endProteinCountPos + 1)) # +1 skip header
			elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
				readingDepth1 = getReadingDepth(GetLine(File, endProteinCountPos + 1)) # +1 skip header
				readingDepth2 = getReadingDepth(GetLine(File, endProteinCountPos + 1 + 1)) # +1 skip header
				readingDepth3 = getReadingDepth(GetLine(File, endProteinCountPos + 2 + 1)) # +1 skip header
				readingDepth = min(readingDepth1, readingDepth2, readingDepth3)
			
		
			if (readingDepth < readMinimalCodonDepth):
				startProteinCountPos = -9999
				break				
			
			startProteinCountPos = endProteinCountPos - (SEGMENTS_DESIRED_LENGTHS[segment] * readPerLine)
			
			if(startProteinCountPos < 0):
				foundStartPositionWithMinimalDepth = True # Accept when M is not in file (because it is the same as having a low readingDepth at this position)
				break
				
			# Check if the expected position has a low readingDepth and accept it.
			# Otherwice, check if it is a M and accept that (can be the case when we found the first M in the consensus, while the correct M is located a few position later)
			if (DNAProteinType == __DATATYPE_AMINOACIDS):
				readingDepth = getReadingDepth(GetLine(File, startProteinCountPos + 1)) # +1 skip header
			elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
				readingDepth1 = getReadingDepth(GetLine(File, startProteinCountPos + 1)) # +1 skip header
				readingDepth2 = getReadingDepth(GetLine(File, startProteinCountPos + 1 + 1)) # +1 skip header
				readingDepth3 = getReadingDepth(GetLine(File, startProteinCountPos + 2 + 1)) # +1 skip header
				readingDepth = min(readingDepth1, readingDepth2, readingDepth3)
			
			if (readingDepth < readMinimalCodonDepth):
				foundStartPositionWithMinimalDepth = True
				break
			else:
				if (DNAProteinType == __DATATYPE_AMINOACIDS):
					if (consensusSequences[startProteinCountPos] == PROTEIN_START_AA.upper()):
						foundStartPositionWithMinimalDepth = True
					else:
						startProteinCountPos = -9999
						break
				elif(DNAProteinType == __DATATYPE_NUCLEOTIDES):
					if(consensusSequences[startProteinCountPos: (startProteinCountPos + readPerLine)] == PROTEIN_START_NT.upper()):
						foundStartPositionWithMinimalDepth = True
					else:
						startProteinCountPos = -9999
						break
			break
			
		else:
			foundStartPositionWithMinimalDepth = True
	
	if(startProteinCountPos == -9999):
		# Read out how many positions in this file equals or exceeds the readMinimalCodonDepth
		perc = GetPercentageAboveProteinReadingDepthThreshold(File, offset)
		writeToWarningLog("File: %s (%s%% above minimal reading depth (%s) used in Variant Detector)" % (File, str(("%.0f" % perc)), rejectOnStopCodonBeforeLine), True)
		return startProteinCountPos
	
	startProteinCountPos = startProteinCountPos + 1
	return startProteinCountPos 
	
# Insert data on a specific position into a line
def insertIntoLine (line, insertIndex, pos):
	return line[:pos] + tabularDivide + insertIndex + line[pos:]

# Insert data on a specific position into a filename
def insertIntoFileName (line, insertIndex, pos):
	return line[:pos] + insertIndex + line[pos:]

# Add count into current file content
def AddProteinOrSegmentCount(File, startProteinCountPos, outputFolder):
	index = (startProteinCountPos - 1) * -1 # Start with negative count until we meet the protein start position
	segment = GetFileSegment(File)

	# Make file to work with, a temp file.
	# So we can save our data into the same named file
	sTempInFileName = ".temp." + File
	
	os.rename(outputFolder + "/" + File, outputFolder + "/" + sTempInFileName)

	if(not(OverWriteOutputFiles) and os.path.exists(outputFolder + "/" + File)):
		QuitAndShowError([("Unable to write to output, since file already exists."),
		("File: " + outputFolder + "/" + File)])
	
	fileOut = open(outputFolder + "/" + File,'w')
	fileOut.truncate() # Empty the file if it already exists
	isHeader = True
	
	# Determine if negative proteinCount lines should be written or not
	DNAProteinType = GetFileDNAProteinType(File)
	deleteNegativeProteinCountedLines = (DNAProteinType == __DATATYPE_AMINOACIDS)
	previousPosition = -1
	for line in open(outputFolder + "/" + sTempInFileName):
		consensus = getConsensus(line)
			 
		writeProteinCountOnPos = find_nth(line,tabularDivide, COLUMNS_ADD_TO_POSITION)
		newLine = ""
		if(isHeader):
			isHeader = False
			newLine = insertIntoLine(line, HEADER_PROTEINCOUNT, writeProteinCountOnPos)
			fileOut.write(newLine)
		else:	
			if (index >= 0 and consensus == stopCodon):
				break

			# check if we have a consecutive numbering.
			# If not; abort.
			# Exception for M2 and NS2 (described in PROTEIN_GAPS dictionary); allow known gap
			currentPosition = getPositionCount(line)
			if(not previousPosition == -1):
				if (currentPosition - previousPosition > 1):
					if(not ((len(PROTEIN_GAPS[segment]) > 0) and (PROTEIN_GAPS[segment][0] == previousPosition and PROTEIN_GAPS[segment][1] == currentPosition))):
						QuitAndShowError([("File does not have consecutive numbering."),
							("Repair or remove file and restart script (note: Provide complete output data from previous step)"),
							("Previous line number: " + str(previousPosition)),
							("Current line number: " + str(currentPosition)),
							("In file: " + File)])
			previousPosition = currentPosition
			
			newLine = insertIntoLine(line, str(index + 1), writeProteinCountOnPos)
			
			if (deleteNegativeProteinCountedLines):
				if (index >= 0):
					fileOut.write(newLine)
			else:
				fileOut.write(newLine)	
				
			index = index + 1
		
	# Clean up
	fileOut.close()
	os.remove(outputFolder + "/" + sTempInFileName)
	
# Add count to all files
def AddProteinAndSegmentCountToFiles(Files, outputFolder):		
	for file in Files:
		Show("Add count to file: " + file)
		consensusSequences = GetConsensusSequence(file)
		startCountPos = GetProteinOrSegmentCountStartPosition(file, consensusSequences)
		if(not(startCountPos == -9999)):
			AddProteinOrSegmentCount(file, startCountPos, outputFolder)
		else:
			os.remove(outputFolder + "/" + file) # Remove this file so we will not see this in the following analysis
			
		
# Initialise arguments
def Initialise():
	args = len(sys.argv)
	if(not (args >= 5)):
		QuitAndShowError([("Not all arguments are provided.")])
		
	global inputFolder
	global inputFolderNT
	global outputFolder
	global outputFolderNT
	global outputFolderAA
	global SEGMENTS_DESIRED_LENGTHS
	global PROTEIN_GAPS
	global _warning_file_handler
	
	inputFolder = sys.argv[1]
	inputFolderNT = sys.argv[2]
	outputFolderNT = sys.argv[3]
	outputFolderAA = sys.argv[4]
	
	outputFolder = outputFolderNT

	
	
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_PB1N40] = __SEGMENTS_DESIRED_LENGTHS_PB1N40
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_PB2] = __SEGMENTS_DESIRED_LENGTHS_PB2
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_PB1] = __SEGMENTS_DESIRED_LENGTHS_PB1
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_PA] = __SEGMENTS_DESIRED_LENGTHS_PA
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_HA] = __SEGMENTS_DESIRED_LENGTHS_HA
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_HAH3] = __SEGMENTS_DESIRED_LENGTHS_HA # Do not use HAH3 since we do not check for the Q, but just for the M
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_NP] = __SEGMENTS_DESIRED_LENGTHS_NP
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_NA] = __SEGMENTS_DESIRED_LENGTHS_NA
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_PB1F2] = __SEGMENTS_DESIRED_LENGTHS_PB1F2
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_PAX] = __SEGMENTS_DESIRED_LENGTHS_PAX
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_MP] = __SEGMENTS_DESIRED_LENGTHS_M1
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_M1] = __SEGMENTS_DESIRED_LENGTHS_M1
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_M2] = __SEGMENTS_DESIRED_LENGTHS_M2
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_NS] = __SEGMENTS_DESIRED_LENGTHS_NS1
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_NS1] = __SEGMENTS_DESIRED_LENGTHS_NS1
	SEGMENTS_DESIRED_LENGTHS[SEGMENTS_NS2] = __SEGMENTS_DESIRED_LENGTHS_NS2

	
	PROTEIN_GAPS[SEGMENTS_PB1N40] = __PROTEIN_GAP_PB1N40
	PROTEIN_GAPS[SEGMENTS_PB2] = __PROTEIN_GAP_PB2
	PROTEIN_GAPS[SEGMENTS_PB1] = __PROTEIN_GAP_PB1
	PROTEIN_GAPS[SEGMENTS_PA] = __PROTEIN_GAP_PA
	PROTEIN_GAPS[SEGMENTS_HA] = __PROTEIN_GAP_HA
	PROTEIN_GAPS[SEGMENTS_HAH3] = __PROTEIN_GAP_HAH3
	PROTEIN_GAPS[SEGMENTS_NP] = __PROTEIN_GAP_NP
	PROTEIN_GAPS[SEGMENTS_NA] = __PROTEIN_GAP_NA
	PROTEIN_GAPS[SEGMENTS_PB1F2] = __PROTEIN_GAP_PB1F2
	PROTEIN_GAPS[SEGMENTS_PAX] = __PROTEIN_GAP_PAX
	PROTEIN_GAPS[SEGMENTS_MP] = __PROTEIN_GAP_MP
	PROTEIN_GAPS[SEGMENTS_M1] = __PROTEIN_GAP_M1
	PROTEIN_GAPS[SEGMENTS_M2] = __PROTEIN_GAP_M2
	PROTEIN_GAPS[SEGMENTS_NS] = __PROTEIN_GAP_NS
	PROTEIN_GAPS[SEGMENTS_NS1] = __PROTEIN_GAP_NS1
	PROTEIN_GAPS[SEGMENTS_NS2] = __PROTEIN_GAP_NS2

	Show("Output path: " + outputFolder)
	
	_warning_file_handler = open(outputFolder + "/" + outputWarningLogFile, 'w')
	_warning_file_handler.truncate() # Empty the file if it already exists

# Return path to dual frame file	
def GetPathToFrameFile(Files, Segment, ProteinType, SecondFrame):
	numberOfFilesFound = 0
	frameNumbersFound = []
	pathToFile = []
	
	# Find out which one is the Fx+1
	for File in Files:
		DNAProteinType = GetFileDNAProteinType(File)
		if (DNAProteinType == ProteinType):
			if(GetFileSegment(File) == Segment):
				frameNumbersFound.append(GetReadingFrameNumber(File))
				pathToFile.append(File)
				numberOfFilesFound = numberOfFilesFound + 1
				
	if (not numberOfFilesFound == 2):
		Show("Segment '" + Segment + "' will be skipped (no dual reading frame found). Found " + str(numberOfFilesFound) + " file(s).")
		return None
		
	# Figure out which one is the second frame
	secondFrameFilePath = ""

	if(SecondFrame):
		# Do we have a f3 - f1 situation
		if(abs(frameNumbersFound[0] - frameNumbersFound[1]) > 1):
			if(frameNumbersFound[0] < frameNumbersFound[1]):
				secondFrameFilePath  = pathToFile[0]
			else:
				secondFrameFilePath  = pathToFile[1]
		elif (frameNumbersFound[0] < frameNumbersFound[1]):
			secondFrameFilePath  = pathToFile[1]
		else:
			secondFrameFilePath  = pathToFile[0]
	else:
		# Do we have a f3 - f1 situation
		if(abs(frameNumbersFound[0] - frameNumbersFound[1]) > 1): 
			if(frameNumbersFound[0] < frameNumbersFound[1]):
				secondFrameFilePath  = pathToFile[1]
			else:
				secondFrameFilePath  = pathToFile[0]
		elif (frameNumbersFound[0] < frameNumbersFound[1]):
			secondFrameFilePath  = pathToFile[0]
		else:
			secondFrameFilePath  = pathToFile[1]
	
		
	return secondFrameFilePath
	
# Change segment inside file content
def AlterSegment(File, add_to_filename, segment, headerName, startColumnCountAtPos, startColumnAAShouldBe, keepOriginalFile):

	# Add headerName to output filename
	Show("Add " + headerName + " count to file: " + File)
	index = 1
	
	# Get output file to work with
	segmentLocationInFileName = File.rfind("_")
	outFile = insertIntoFileName(File, add_to_filename, segmentLocationInFileName)
	
	if(not(OverWriteOutputFiles) and os.path.exists(outputFolder + "/" + outFile)):
		QuitAndShowError([("Unable to write to output, since file already exists."),
		("File: " + outputFolder + "/" + File)])
	
	fileOut = open(outputFolder + "/" + outFile, 'w')
	fileOut.truncate() # Empty the file if it already exists
	isHeader = True
	
	# Determine if negative currentPos lines should be written or not
	proteinStarted = False
	proteinStartedAtPos = -1
	for line in open(inputFolder + "/" + File):
		writeColumnOnPos = find_nth(line, tabularDivide, COLUMNS_ADD_TO_POSITION)
		
		if(isHeader):
			isHeader = False
			newLine = insertIntoLine(line, headerName, writeColumnOnPos)
			fileOut.write(newLine)
		else:
			currentPos = getPositionCount(line)
			consensus = getConsensus(line)
			
			readingDepth = getReadingDepth(line)

			if(not(proteinStarted) and consensus == PROTEIN_START_AA and (readingDepth >= readMinimalCodonDepth)):
				proteinStarted = True
				proteinStartedAtPos = currentPos - 1

			newLine = ""
			if (proteinStarted == True and (currentPos - proteinStartedAtPos)>= startColumnCountAtPos):
				newLine = insertIntoLine(line, str(index), writeColumnOnPos)
				
				# Check if first consensus is as expected
				if(index == 1):
					consensus = getConsensus(line)
					if(not consensus == startColumnAAShouldBe):
						Show("WARNING: " + segment + ": " + str(startColumnAAShouldBe) + " not at " + str(startColumnCountAtPos) + ", " + segment + " numbering incorrect. Found consensus: " + consensus + " on " + str(currentPos))
				index = index + 1
			else:
				newLine = insertIntoLine(line, "", writeColumnOnPos)
			
			fileOut.write(newLine)
			
	fileOut.close()
	
	# Keep original file if required
	if(keepOriginalFile):
		copyfile(inputFolder + "/" + File, outputFolder + "/" + File)
			
# Change specific segment files
def AlterSegmentFiles(Files):

	# HA
	Show("Alter HA segment files...")
	for File in Files:
		DNAProteinType = GetFileDNAProteinType(File)
		if (DNAProteinType == __DATATYPE_AMINOACIDS):
			if(GetFileSegment(File) == SEGMENTS_HA):
				AlterSegment(File, SEGMENTS_HA_H3_ADD_TO_FILENAME, SEGMENTS_HA, HEADER_H3, SEGMENT_HA_START_H3_COUNT_AT_POS, SEGMENT_HA_START_H3_AA_SHOULD, False)
			
	# PB1
	Show("Alter PB1 segment files...")
	firstFrameFilePath = GetPathToFrameFile(Files, SEGMENTS_PB1, __DATATYPE_AMINOACIDS, False)
	if(firstFrameFilePath is not None):
		AlterSegment(firstFrameFilePath, SEGMENTS_PB1_N40_ADD_TO_FILENAME, SEGMENTS_PB1, HEADER_N40, SEGMENT_PB1_START_N40_COUNT_AT_POS, SEGMENT_PB1_START_N40_AA_SHOULD, True)
	
# Check for lengths of protein files
def CheckAllProteinFiles(proteinFiles):
	
	Show ("Check if all proteins are as expected")
	for file in proteinFiles:
		Show("Check file: " + file)
		consensus = GetConsensusSequence(file)
		segment = GetFileSegment(file)

		checkSegmentsForDesiredLength = [	SEGMENTS_PB2, 
											SEGMENTS_PB1, 
											SEGMENTS_PA, 
											SEGMENTS_HA, 
											SEGMENTS_PB1N40, 
											SEGMENTS_NP, 
											SEGMENTS_NA, 
											SEGMENTS_PB1F2, 
											SEGMENTS_PAX,
											SEGMENTS_M1,
											SEGMENTS_M2,
											SEGMENTS_NS1,
											SEGMENTS_NS2]
											
		if(segment in checkSegmentsForDesiredLength): 
			if(not(len(consensus) == SEGMENTS_DESIRED_LENGTHS[segment])):
				writeToWarningLog("%s: Length (%d) not as expected: %d" % (segment, len(consensus), SEGMENTS_DESIRED_LENGTHS[segment]))
				
		if(segment == SEGMENTS_PB1N40):
			pos = getPositionCount(GetLine(file, 1))
			if(not(pos == SEGMENTS_PB1_MN40_POS)):
				writeToWarningLog("PB1-N40: M at unexpected position")

		if(segment == SEGMENTS_HA):				
			c = getConsensus(GetLine(file, 17))
			if(not(c == SEGMENT_HA_START_H3_AA_SHOULD)):
				writeToWarningLog("HA: Q not at 17, HA numbering incorrect")

		if(segment == SEGMENTS_PB1F2):
			pos = getPositionCount(GetLine(file, 1))
			if(not(pos == SEGMENTS_PB1_SECOND_M_POS)):
				writeToWarningLog("PB1-F2: M at unexpected position")

		if(segment == SEGMENTS_PAX):
			pos = getPositionCount(GetLine(file, 1))
			if(not(pos == SEGMENTS_PA_M_POS_FRAME_FIRST)):
				writeToWarningLog("PA-X: M at unexpected position")

		if(segment == SEGMENTS_M2):
			pos = getPositionCount(GetLine(file, 1))
			if(not(pos == SEGMENTS_M1_START_PROTEINCOUNT_POS_FRAME1)):
				writeToWarningLog("M2: M at unexpected position")

		if(segment == SEGMENTS_NS2):
			pos = getPositionCount(GetLine(file, 1))
			if(not(pos == SEGMENTS_NS1_START_PROTEINCOUNT_POS_FRAME1)):
				writeToWarningLog("NS2: M at unexpected position")
				
		if(segment == SEGMENTS_PAX):
			consensus = getConsensus(GetLine(file, SEGMENTS_PA_V_POS_FRAME_SECOND - getPositionCount(GetLine(file, 1)) + 1))
			if(not(consensus == "V")):
				writeToWarningLog(("PA-X: V not at expected position: %d" % (SEGMENTS_PA_V_POS_FRAME_SECOND)))
				
	for line in open(outputFolder + "/" + outputWarningLogFile):
		Show(line)
	
# Concat two files and remove the second source file
def CreateCombinedProteinFileAndRemoveSecondFile(Files, Segment, SegmentNameFirst, SegmentNameSecond, FirstFrameStart, FirstFrameEnd, SecondFrameStart):
	
	firstFrameFilePath = GetPathToFrameFile(Files, Segment, __DATATYPE_AMINOACIDS, False)
	if(firstFrameFilePath is None):
		return
		
	secondFrameFilePath = GetPathToFrameFile(Files, Segment, __DATATYPE_AMINOACIDS, True)
	if(secondFrameFilePath is None):
		return

	Show("Using the following reading frame files:")
	Show("-" +  firstFrameFilePath)
	Show("-" + secondFrameFilePath)
		
		
	outFile = firstFrameFilePath.replace(Segment, SegmentNameSecond)
	if(not(OverWriteOutputFiles) and os.path.exists(outputFolder + "/" + outFile)):
		QuitAndShowError([("Unable to write to output, since file already exists."),
		("File: " + outputFolder + "/" + File)])
	
	fileOut = open(outputFolder + "/" + outFile, 'w')
	fileOut.truncate() # Empty the file if it already exists

	# Get part from first file
	isHeader = True
	for line in open(inputFolder + "/" + firstFrameFilePath):
		if(isHeader):
			fileOut.write(line)
			isHeader = False
		else:
			proteinCount = getPositionCount(line)
			
			if (proteinCount >= FirstFrameStart):
				if(proteinCount <= FirstFrameEnd):
					fileOut.write(line)
				else:
					break
				
	# Get part from second file
	isHeader = True

	for line in open(inputFolder + "/" + secondFrameFilePath):
		if(isHeader):
			isHeader = False
		else:
			proteinCount = getPositionCount(line)
			if(proteinCount >= SecondFrameStart):
			
				consensus = getConsensus(line)
				if(consensus == stopCodon):
					break
				
				fileOut.write(line)
				
	# Copy first segment file
	newFileName_SegmentNameFirst = firstFrameFilePath.replace(Segment, SegmentNameFirst)
	copyfile(inputFolder + "/" + firstFrameFilePath, outputFolder + "/" + newFileName_SegmentNameFirst)
	
# Rename file PB1 to PB1F2 
def ChangePB1F2File(Files):
	secondFrameFilePath = GetPathToFrameFile(Files, SEGMENTS_PB1, __DATATYPE_AMINOACIDS, True)
	if(secondFrameFilePath is not None):
		newFileName = secondFrameFilePath.replace("PB1", "PB1F2")
		copyfile(inputFolder + "/" + secondFrameFilePath, outputFolder + "/" + newFileName)
	
# Change listed file names
def ChangeFileNames(Files):
	Show ("Rename PB1 frame2 file to PB1F2")
	ChangePB1F2File(Files)

# Copy non altered segment file
def CopyNonAlteredFile(Files, segment):
	for File in Files:
		DNAProteinType = GetFileDNAProteinType(File)
		if (DNAProteinType == __DATATYPE_AMINOACIDS):
			if(GetFileSegment(File) == segment):
				Show("Copy file to output folder: " + File)
				copyfile(inputFolder + "/" + File, outputFolder + "/" + File)
				
# Copy listed non altered file
def CopyNonAlteredFiles(Files):
	CopyNonAlteredFile(Files, SEGMENTS_NP)
	CopyNonAlteredFile(Files, SEGMENTS_NA)
	CopyNonAlteredFile(Files, SEGMENTS_PB2)

# Build new protein files, each consist of 2 source segment files
def CreateNewProteinFiles(Files):
	Show ("Create NS1 and NS2 file")
	CreateCombinedProteinFileAndRemoveSecondFile(Files, SEGMENTS_NS, "NS1", "NS2", SEGMENTS_NS1_START_PROTEINCOUNT_POS_FRAME1, SEGMENTS_NS2_STOP_PROTEINCOUNT_POS_FRAME1, SEGMENTS_NS2_START_PROTEINCOUNT_POS_FRAME2)
	Show ("Create M1 and M2 file")
	CreateCombinedProteinFileAndRemoveSecondFile(Files, SEGMENTS_MP, "M1", "M2",SEGMENTS_M1_START_PROTEINCOUNT_POS_FRAME1,  SEGMENTS_M2_STOP_PROTEINCOUNT_POS_FRAME1, SEGMENTS_M2_START_PROTEINCOUNT_POS_FRAME2)
	Show ("Create PA and PAX file")
	CreateCombinedProteinFileAndRemoveSecondFile(Files, SEGMENTS_PA, "PA", "PAX", SEGMENTS_PA_START_PROTEINCOUNT_POS_FRAME1, SEGMENTS_PA_STOP_PROTEINCOUNT_POS_FRAME1, SEGMENTS_PA_START_PROTEINCOUNT_POS_FRAME2)
	
# Removing any file which isn't AA typed
def RemoveNonSameTypedFiles(Files, remainFileType):
	remainingFiles = []
	
	for file in Files:
		if (GetFileDNAProteinType(file) == remainFileType):
			remainingFiles.append(file)
	
	Debug("Found %s typed files: %s" % (remainFileType, str(len(remainingFiles))))
	return remainingFiles

# Get first frame file
def SelectFirstFrameFile(file, similarFiles, dataType):
	if (len(similarFiles) == 0):
		return file
		
	if (len(similarFiles) > 1):
		QuitAndShowError([("UNEXPECTED BEHAVIOUR. SelectFirstFrameFile. Expect only 2 readingframes per segment."),
		("Number of reading frames provided: " + str(len(similarFiles) + 1)),
		("One of the readingframe files: " + file)])
		
	files = []
	files.append(file)
	files.append(similarFiles[0])
	firstFrameFilePath = GetPathToFrameFile(files, GetFileSegment(file), dataType, False)
	return firstFrameFilePath

	
# Get first frame files
def SelectFirstFrameFiles(Files, dataType):
	ProcessedFiles = [] 
	copiedFiles = []
	for file in Files:
		if (file in ProcessedFiles):
			continue
		# don't use this file again
		ProcessedFiles.append(file)
		
		# Find similar files
		similarFiles = FindSimilarFiles(Files, file)
		ProcessedFiles.extend(similarFiles)
		
		# Copy the correct reading frame file or place a warning file
		fileToAppend = SelectFirstFrameFile(file, similarFiles, dataType)
		if (fileToAppend is not None):
			copiedFiles.append(SelectFirstFrameFile(file, similarFiles, dataType))
	
	return copiedFiles


# Main method
def Main():
	global inputFolder
	global outputFolder
	global COLUMNS_CONSENSUS
	
	ShowSplash()
	Initialise()
	
	Files = GetFileDirectPaths(inputFolder)
	
	Show("")
	Show("Processing amino acids files...")
	proteinFiles = RemoveNonSameTypedFiles(Files, __DATATYPE_AMINOACIDS)
	
	FilesNT = GetFileDirectPaths(inputFolderNT)
	nucleotideFiles = RemoveNonSameTypedFiles(FilesNT, __DATATYPE_NUCLEOTIDES)
	
	if(len(proteinFiles) + len(nucleotideFiles) == 0):
		QuitAndShowError([("No files to work with.")])
	
	if(len(proteinFiles)  > 0):
		Show("Alter specific segment files...")
		# Add meta data to specific segment files
		# - Creates new HA file in output
		AlterSegmentFiles(proteinFiles)
		
		# Show warnings if specific rules are not met for specific segments
		# CheckSegmentFiles(proteinFiles)
		CreateNewProteinFiles(proteinFiles)
		ChangeFileNames(proteinFiles)
		CopyNonAlteredFiles(proteinFiles)
	
	#################
	## AminoAcids  ##
	#################
	Show("")
	Show("Add protein count to all files...")
	inputFolderOld = inputFolder
	inputFolder = outputFolder
	FilesToAddCount = GetFileDirectPaths(inputFolder)

	writeToWarningLog("WARNING. These AA files do not contain a consensus sequence which describes the start of the protein and the stopcodon location with enough depth (%s). These file will be ignored:" % (str(readMinimalCodonDepth)))
	AddProteinAndSegmentCountToFiles(FilesToAddCount, outputFolder)
	
	# Check end files
	filesToCheck = GetFileDirectPaths(outputFolder)
	CheckAllProteinFiles(filesToCheck)
	
	#################
	## Nucleotides ##
	#################
	Show("")
	Show("Add segment count to all files...")	

	inputFolder = inputFolderNT
	outputFolder = outputFolderAA

	NTFirstFrameFilesOnly = SelectFirstFrameFiles(nucleotideFiles, __DATATYPE_NUCLEOTIDES)
	if(len(NTFirstFrameFilesOnly)  > 0):
		Show("Copy " + str(len(NTFirstFrameFilesOnly)) + " nucleotide files to output...")
		CopyFilesToOutput(inputFolder, NTFirstFrameFilesOnly, outputFolder)
		
	#inputFolder = outputFolderAA
	FilesToAddCount = GetFileDirectPaths(outputFolder)
	writeToWarningLog("WARNING. These nt files do not contain a consensus sequence which describes the start of the protein and the stopcodon location with enough depth (%s). These file will be ignored:" % (str(readMinimalCodonDepth)))
	AddProteinAndSegmentCountToFiles(FilesToAddCount, outputFolder)

	_warning_file_handler.close()
	Show("Program ready!")

Main()