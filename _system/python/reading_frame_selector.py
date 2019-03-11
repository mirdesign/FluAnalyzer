#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# Reading frame selector                                    #
#############################################################
from __future__ import print_function
import sys 
import os
import glob
from shutil import copyfile

# Load common shared variables and methods
execfile(os.path.join(os.path.dirname(__file__),"_global_constants.py"))

inputFolder = ""
inputFolderAA = ""
inputFolderNT = ""
outputFolder = ""
outputFolderNucleotides = ""


# Show some startup information
def ShowSplash():
	Debug ("")
	Debug ("--------------------------------")
	Debug ("Tool: Readingframe selector")
	Debug ("Purpose: Parse Multiple influenza reading frame files and picks the correct rading frame for each segment")
	Debug ("version: 1.0")
	Debug ("Contact: m.pater@amc.nl | mpater@mirdesign.nl - Maarten Pater")
	Debug ("--------------------------------")
	Debug ("Usage: python rfs.py folder_with_reading_frame_files output_folder_proteins output_folder_nucleotides [readMinimalCodonDepth = 100] [rejectOnStopCodonBeforeLine = 100]")
	Debug ("Expected Filename strategy: [CUSTOM]_[CUSTOM]_T[0..n]_[LOCATION]_f[READINGFRAME]-[segment]_[" + __DATATYPE_AMINOACIDS + "/" + __DATATYPE_NUCLEOTIDES + "].ext")
	Debug ("--------------------------------")


# Get AA or nt type information from filename
def GetFileDNAProteinType(file):
	return file[file.rfind(FILENAMES_DATATYPE_LASTDIVIDE_CHAR) + len(FILENAMES_DATATYPE_LASTDIVIDE_CHAR): file.rfind(FILENAMES_DATATYPE_LASTDIVIDE_CHAR) + len(__DATATYPE_AMINOACIDS) + len(FILENAMES_DATATYPE_LASTDIVIDE_CHAR)]

# Removing any file from the list which isn't AA-typed
def RemoveNonSameTypedFiles(Files, remainFileType):
	remainingFiles = []
	
	for file in Files:
		if (GetFileDNAProteinType(file) == remainFileType):
			remainingFiles.append(file)
	
	Debug("Found %s typed files: %s" % (remainFileType, str(len(remainingFiles))))
	return remainingFiles
	
# Find inside files list similar (based on segment) files based on file
def FindSimilarFiles(Files, file):
	similarFiles = []
	
	currentSegment = GetFileSegment(file)
	Debug("Analysing segment: %s" % currentSegment)
	for similarFile in Files:
		if (similarFile == file):
			continue	
		
		if(GetFileSegment(similarFile) == currentSegment):
			similarFiles.append(similarFile)
		
	return similarFiles
	
# Check if the sequence is in a valid reading frame
def HasCorrectReadingFrame(File):

	Show("Check if file has correct readingframe: %s" % File, False)
	bStartCodonFound = False
	cntPositionCheckForStopCodon = 0
	cntStopCodonsFound = 0
	
	
	for line in open(inputFolder + "/" + File):
		startConsensus = find_nth(line,tabularDivide, 2)
		endConsensus = find_nth(line,tabularDivide, 3)
		consensus = line[startConsensus + 1: endConsensus]
		
		# Start looking for stop codons only after start codon has been found
		if ((not bStartCodonFound) and (consensus.upper() == PROTEIN_START_AA.upper())):
			bStartCodonFound = True
		
		if (bStartCodonFound):

			if (consensus == stopCodon):
				cntStopCodonsFound = cntStopCodonsFound + 1
				
				# Stop on exceeding maximum number of codons per file
				if (cntStopCodonsFound == rejectOnStopCodonNumberFound):
					Show (" [No]")
					return False
				
				if (cntPositionCheckForStopCodon <= rejectOnStopCodonBeforeLine):	
					# Check if this read is high enough to trust
					readingDepth = getReadingDepth(line)
					
					# Stop on stop codon found before threshold position with a high enough reading depth
					if(readingDepth >= readMinimalCodonDepth):
						Show (" [No]")
						return False #We found a stopCodon with a high readingDepth, so consider this a "true" read

			cntPositionCheckForStopCodon = cntPositionCheckForStopCodon + 1
	Show (" [Yes]")
	return True	
	

# Analyse files and determine if it they contain a sequence with a valid reading frame 
def GetFilesWithCorrectReadingFrame(file, similarFiles):
	retVal = []
	isCorrectReadingFrame = []
	pathToFile = []
	readingFrameID = []
	segment = GetFileSegment(file)
	lastCorrectSecondReadingFrameID = -1
	lastCorrectReadingFrameID = -1
	
	# First file
	pathToFile.append(file)
	isCorrectReadingFrame.append(HasCorrectReadingFrame(file))
	Show(file)
	readingFrameID.append(GetReadingFrameNumber(file))
	
	
	for sFile in similarFiles:
		pathToFile.append(sFile)
		isCorrectReadingFrame.append(HasCorrectReadingFrame(sFile))
		Show(file)
		readingFrameID.append(GetReadingFrameNumber(sFile))
		

	fileCounter = 0	
	for crf in isCorrectReadingFrame:
		if(crf):
			retVal.append(pathToFile[fileCounter])
			lastCorrectReadingFrameID = readingFrameID[fileCounter]
			
			if(	segment == SEGMENTS_PB1 or 
				segment == SEGMENTS_PA or
				segment == SEGMENTS_NS or
				segment == SEGMENTS_MP):
				
				# Find next readingframe file
				lastCorrectSecondReadingFrameID = lastCorrectReadingFrameID + 1
				if (lastCorrectSecondReadingFrameID > 3): # 3 codons; so 3 readingframes possible
					lastCorrectSecondReadingFrameID = 1
				
				rfID_FileCounter = 0
				for rfID in readingFrameID:
					if (rfID == lastCorrectSecondReadingFrameID):
						retVal.append(pathToFile[rfID_FileCounter])
						break
				
					rfID_FileCounter = rfID_FileCounter + 1
			break
		fileCounter = fileCounter + 1
		
	if (lastCorrectReadingFrameID == -1):
		Show ("lastCorrectReadingFrameID: Not found")
	else:
		Show ("lastCorrectReadingFrameID: " + str(lastCorrectReadingFrameID))
	
	if (lastCorrectSecondReadingFrameID == -1):
		Show ("lastCorrectSecondReadingFrameID: Not found")
	else:
		Show ("lastCorrectSecondReadingFrameID: " + str(lastCorrectSecondReadingFrameID))

	# Check if second correct reading frame is indeed the upfollowing reading frame.
	# Only valid compbination: f1 f2, f2 f3, f3 f1
	if(not (lastCorrectReadingFrameID == -1 or lastCorrectSecondReadingFrameID == -1)):
		if(not(
		(lastCorrectReadingFrameID == 1 and lastCorrectSecondReadingFrameID == 2) or 
		(lastCorrectReadingFrameID == 2 and lastCorrectSecondReadingFrameID == 3) or 
		(lastCorrectReadingFrameID == 3 and lastCorrectSecondReadingFrameID == 1))):
			Show("WARNING! Found correct reading frame to be F" + str(lastCorrectReadingFrameID) + " and correct Second reading frame to be F" + str(lastCorrectSecondReadingFrameID) + ". This is illegal.")
			Show("WARNING! Since this is just a WARNING, we'll continue using these files.")
			#WARNING! Found correct reading frame to be F1 and correct Second reading frame to be F3. This is illegal.

	return retVal

	
# Write a log file which indicates that no file could be found with a correct reading frame
def WriteNoCorrectReadingFrame(file):
	global outputFolder
	open(outputFolder + "/" + file + "__NoCorrectReadingFrameFound_" + GetFileSegment(file) + ".log", 'a')
	
# Write a log file which indicates that multiple files were found with a correct reading frame
def WriteUnableToPickCorrectReadingFrame(file):
	global outputFolder
	open(outputFolder + "/" + file + "__MultipleCorrectReadingFrames_" + GetFileSegment(file) + ".log", 'a')

# Write a log file which indicates that multiple files were found with a correct reading frame
def WriteMissingSecondReadingFrame(file):
	global outputFolder
	open(outputFolder + "/" + file + "__WriteMissingSecondReadingFrame_" + GetFileSegment(file) + ".log", 'a')

	
# Determine if one or multiple files can be found with a correct reading frame
def DetermineCorrectReadingFrame(file, similarFiles):
	filesWithCorrectReadingFrame = GetFilesWithCorrectReadingFrame(file, similarFiles)
	segment = GetFileSegment(file)
	Show("Correct reading frame #files: %s" % str(len(filesWithCorrectReadingFrame)))
	Show("")
	
	copiedFiles = []
	
	if(	segment == SEGMENTS_PB1 or
		segment == SEGMENTS_PA or
		segment == SEGMENTS_NS or
		segment == SEGMENTS_MP):
	
		if (len(filesWithCorrectReadingFrame) == 0):
			WriteNoCorrectReadingFrame(file)
		elif (len(filesWithCorrectReadingFrame) == 1):
			for f in filesWithCorrectReadingFrame :
				WriteMissingSecondReadingFrame(f)
		elif (len(filesWithCorrectReadingFrame) == 2):
			for f in filesWithCorrectReadingFrame :
				CopyFileToOutput(inputFolder, f, outputFolder)
				copiedFiles.append(f)
		else:
			for f in filesWithCorrectReadingFrame :
				WriteUnableToPickCorrectReadingFrame(f)
	else:
		if (len(filesWithCorrectReadingFrame) == 0):
			WriteNoCorrectReadingFrame(file)
		elif (len(filesWithCorrectReadingFrame) == 1):
			CopyFileToOutput(inputFolder, filesWithCorrectReadingFrame[0], outputFolder)
			copiedFiles.append(filesWithCorrectReadingFrame[0])
		else:
			WriteUnableToPickCorrectReadingFrame(file)
	
	return copiedFiles
	
# Process found files
def ProcessFiles(Files):
	if(len(Files) == 0):
		Show ("No files to process. Quit program")
		sys.exit(0)
		
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
		copiedFiles.extend(DetermineCorrectReadingFrame(file, similarFiles))
	return copiedFiles

# Copy nucleotide files with matching readingframes per segment
def copyNucleotideFiles(processedProteinFiles):
	for proteinFile in processedProteinFiles:
		# Convert into Nucleotide file
		Show ("Search for nucleotide pair file for: " + proteinFile)
		postType = "_"
		nucleotideFile = proteinFile.replace(postType + __DATATYPE_AMINOACIDS, postType + __DATATYPE_NUCLEOTIDES)
		
		# Copy nucleotide file to output
		Show("Copy nucleotide pair file to output: " + nucleotideFile)
		CopyFileToOutput(inputFolder, nucleotideFile, outputFolderNucleotides)

# Initialise arguments
def Initialise():
	args = len(sys.argv)
	if(not (args >= 5)):
		Debug("ERROR Not all arguments are provided")
		sys.exit(1)
		
	global inputFolderAA
	global inputFolderNT
	global outputFolder
	global outputFolderNucleotides
	global readMinimalCodonDepth
	global rejectOnStopCodonBeforeLine
	
	inputFolderAA = sys.argv[1]
	inputFolderNT = sys.argv[2]
	outputFolder = sys.argv[3]
	outputFolderNucleotides = sys.argv[4]
		
	if(args > 5):
		readMinimalCodonDepth = int(sys.argv[5])
	if(args > 6):
		rejectOnStopCodonBeforeLine = int(sys.argv[6])
	
# Main method
def Main():
	global inputFolder
	
	ShowSplash()
	Initialise()
	
	inputFolder = inputFolderAA
	Files = GetFileDirectPaths(inputFolder)
	
	### AA ###
	# Get list of aminoacid typed files
	Show("")
	Show("Processing amino acids files...")
	proteinFiles = RemoveNonSameTypedFiles(Files, __DATATYPE_AMINOACIDS)
	# Process protein files
	processedProteinFiles = ProcessFiles(proteinFiles)
	
	### NT ###
	# Copy nucleotide files
	Show("Copy nucleotide files..")
	inputFolder = inputFolderNT
	
	copyNucleotideFiles(processedProteinFiles)
	
	Show("")
	Show ("Program ready.")
Main()