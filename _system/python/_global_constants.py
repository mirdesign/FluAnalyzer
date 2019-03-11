#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
#  All global used methods + variables                      #
#############################################################
# Maarten Pater (m.pater@amc.nl | mpater@mirdesign.nl)
from __future__ import print_function
import sys 
import os



# Constants Global
debug = True
deepDebug = True
showMessages = True
OverWriteOutputFiles = True
__DATATYPE_AMINOACIDS = "AA" # Can only be 2 chars!
__DATATYPE_NUCLEOTIDES = "NA" # Can only be 2 chars!
__FILENAMES_READINGFRAME_START = "_f"
FILENAMES_DATATYPE_LASTDIVIDE_CHAR = "_"
PROTEIN_START_AA = "M"
PROTEIN_START_NT = "ATG"
tabularDivide =  "\t"


stopCodon = "*"
stopCodon_nt_1 = "TAA"
stopCodon_nt_2 = "TAG"
stopCodon_nt_3 = "TGA"

readMinimalCodonDepth = 10
rejectOnStopCodonBeforeLine = 100
rejectOnStopCodonNumberFound = 7

SEGMENTS_HA = "HA"
SEGMENTS_HAH3 = "HAH3"
SEGMENTS_NP = "NP"
SEGMENTS_NA = "NA"
SEGMENTS_PA = "PA"
SEGMENTS_PAX = "PAX"
SEGMENTS_PB1 = "PB1"
SEGMENTS_PB1F2 = "PB1F2"
SEGMENTS_PB1N40 = "PB1N40"
SEGMENTS_PB2 = "PB2"
SEGMENTS_MP = "MP"
SEGMENTS_M1 = "M1"
SEGMENTS_M2 = "M2"
SEGMENTS_NS = "NS"
SEGMENTS_NS1 = "NS1"
SEGMENTS_NS2 = "NS2"


__SEGMENTS_DESIRED_LENGTHS_PB1N40 = 718
__SEGMENTS_DESIRED_LENGTHS_PB2 = 759
__SEGMENTS_DESIRED_LENGTHS_PB1 = 757
__SEGMENTS_DESIRED_LENGTHS_PA = 716
__SEGMENTS_DESIRED_LENGTHS_HA = 566
__SEGMENTS_DESIRED_LENGTHS_HAH3 = 566 - 16
__SEGMENTS_DESIRED_LENGTHS_NP = 498
__SEGMENTS_DESIRED_LENGTHS_NA = 469
#__SEGMENTS_DESIRED_LENGTHS_MP = not set, same as M1
__SEGMENTS_DESIRED_LENGTHS_PB1F2 = 90
__SEGMENTS_DESIRED_LENGTHS_PAX = 252
__SEGMENTS_DESIRED_LENGTHS_M1 = 252
__SEGMENTS_DESIRED_LENGTHS_M2 = 97
__SEGMENTS_DESIRED_LENGTHS_NS1 = 230
__SEGMENTS_DESIRED_LENGTHS_NS2 = 121


# Only 1 gap per protein is allowed.
# Describe gap as followed: [last position before gap, first position after gap]
# i.e.: If your protein looks like 1234-89, the gap would be described as [4, 8]
__PROTEIN_GAP_PB1N40 = []
__PROTEIN_GAP_PB2 = []
__PROTEIN_GAP_PB1 = []
__PROTEIN_GAP_PA = []
__PROTEIN_GAP_HA = []
__PROTEIN_GAP_HAH3 = []
__PROTEIN_GAP_NP = []
__PROTEIN_GAP_NA = []
__PROTEIN_GAP_PB1F2 = []
__PROTEIN_GAP_PAX = []
__PROTEIN_GAP_MP = []
__PROTEIN_GAP_M1 = []
__PROTEIN_GAP_M2 = [17, 247]
__PROTEIN_GAP_NS = []
__PROTEIN_GAP_NS1 = []
__PROTEIN_GAP_NS2 = [17, 176]





SEGMENTS_HA_H3_ADD_TO_FILENAME = "H3"

SEGMENT_HA_START_H3_COUNT_AT_POS = 17 # File Count
SEGMENT_HA_START_H3_AA_SHOULD = "Q"

SEGMENT_PB1_START_N40_COUNT_AT_POS = 40 # protein Count (since method knows where M starts)- was 31
SEGMENT_PB1_START_N40_AA_SHOULD = "M"
SEGMENTS_PB1_MN40_POS = 48 # File Count
SEGMENTS_PB1_SECOND_M_POS = 40 # File Count
SEGMENTS_PB1F2_SECOND_M_POS = 31 # File Count
SEGMENTS_PB1_STOPCODON_POS = 178 # File Count
SEGMENTS_PB1_N40_ADD_TO_FILENAME = "N40"	

SEGMENTS_PA_M_POS_FRAME_FIRST = 9
SEGMENTS_PA_V_POS_FRAME_SECOND = 200 
SEGMENTS_PA_STOPCODON_1_POS = 41 + 200 # 200 is from 1st reading frame
SEGMENTS_PA_STOPCODON_2_POS = 61 + 200 # 200 is from 1st reading frame

SEGMENTS_PA_START_PROTEINCOUNT_POS_FRAME1 = 9 # File count
SEGMENTS_PA_STOP_PROTEINCOUNT_POS_FRAME1 = 199 # File count
SEGMENTS_PA_START_PROTEINCOUNT_POS_FRAME2 = 200 # File count

SEGMENTS_NS1_START_PROTEINCOUNT_POS_FRAME1 = 9 # File count
SEGMENTS_NS2_STOP_PROTEINCOUNT_POS_FRAME1 = 17 # File count
SEGMENTS_NS2_START_PROTEINCOUNT_POS_FRAME2 = 176 # File count

SEGMENTS_M1_START_PROTEINCOUNT_POS_FRAME1 = 9 # File count
SEGMENTS_M2_STOP_PROTEINCOUNT_POS_FRAME1 = 17 # File count
SEGMENTS_M2_START_PROTEINCOUNT_POS_FRAME2 = 247 # File count
	
COLUMNS_CONSENSUS = 2
COLUMNS_READINGDEPTH = 3
__OUTPUT_ROW_SEPERATOR = ";"
__OUTPUT_NOT_AVAILABLE_VALUE = "NA"
__INPUT_ROW_SEPERATOR = ";"
__FILE_VARIANT_OVERVIEW_FILE_EXTENTION = "csv"
		

# Get reading frame number from filename
def GetReadingFrameNumber(file):
	return int(file[file.find(__FILENAMES_READINGFRAME_START) + len(__FILENAMES_READINGFRAME_START): file.find("-")])
	
def ShowError(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    Show(*args)
	
# Show error message and quit program
def QuitAndShowError(ErrorMessages):

	headerMessage = "  ERROR  "
	prefixMessage = "ERROR "
	
	# Calculate length of header + footer
	longestMessage = 0
	for message in ErrorMessages:
		if (len(message) + len("ERROR ") > longestMessage):
			longestMessage = len(prefixMessage) + len(message)
		
	numAst = (int)(round(longestMessage - len (headerMessage) + 0.5) / 2)
	
	# Header
	ShowError("\n" + ("*" * numAst) + headerMessage + ("*" * numAst) + "\n")
	
	# Show error message
	for message in ErrorMessages:
		ShowError(prefixMessage + message)
	
	# Footer
	ShowError("\n" + ("*" * numAst) + headerMessage + ("*" * numAst) + "\n")
	sys.exit(1)	
	
# Show information
def Show(message, newLine = True):
	if showMessages:
		if(newLine):
			print (message)
		else:
			print (message,)
		
# Show debug information
def Debug(DebugInfo):
    if debug:
        print (DebugInfo)

# Show deep debug information
def DeepDebug(DebugInfo):
    if deepDebug:
        print (DebugInfo)



# Find a string, n-th occurence
def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start	
	
def getReadingDepth(line):
	startReadingDepth = find_nth(line,tabularDivide, COLUMNS_READINGDEPTH)
	endReadingDepth = find_nth(line,tabularDivide, COLUMNS_READINGDEPTH + 1)
	return int(line[startReadingDepth + 1: endReadingDepth])
	
# Get a list of all files inside the provided folder
def GetFileDirectPaths (Folder):
	Debug("Reading files from: %s " % Folder)
	if (not os.path.isdir(Folder)): 
		Debug("Provided path does not exist: %s " % Folder)
		sys.exit(1)
		
	retval = []
	files = []
	files.extend(os.listdir(Folder))
	
	# Only get files
	for f in files:
		if (os.path.isfile(Folder + "/" + f)):
			if(not(f[-3:] == "log")):
				retval.append(f)
			
	Debug("Found files to process: %s" % str(len(retval)))
	return retval
	
# Get segment name of influenza
def GetFileSegment(file):
	global FILENAMES_DATATYPE_LASTDIVIDE_CHAR
	beforeDNAProteinTypeInfo = file.rfind(FILENAMES_DATATYPE_LASTDIVIDE_CHAR)
	strippedDNAProteinFileName = file[:beforeDNAProteinTypeInfo]
	return strippedDNAProteinFileName[strippedDNAProteinFileName.rfind("-") + 1:]
	
    		
# Get data of a specific column from a tabular divided file
def GetDataFromFileContent(line, COL, divide = None):
	if divide is None:
		divide = __OUTPUT_ROW_SEPERATOR # Default value

	startData = 0
	startDataSkipDividerSize = 0
	if (COL > 0):
		startData = find_nth(line, divide, COL)
		startDataSkipDividerSize = 1
		
	endData = find_nth(line, divide, COL + 1)
	return line[startData + startDataSkipDividerSize: endData]
				
# Find inside files list similar files based on file
def FindSimilarFiles(Files, file):
	similarFiles = []
	
	currentSegment = GetFileSegment(file)
	for similarFile in Files:
		if (similarFile == file):
			continue	
		
		if(GetFileSegment(similarFile) == currentSegment):
			similarFiles.append(similarFile)
		
	return similarFiles
    
# Copy file to output folder
def CopyFileToOutput(inputFolder, fileName, outputfolder = ""):

	srcFile = os.path.join(inputFolder , fileName)
	outFile = os.path.join(outputfolder, fileName)

	if(not(os.path.isfile(srcFile))):
		QuitAndShowError([("CopyFileToOutput: File does not exist."),
							("File provided: " + srcFile)]) # ERROR

	copyfile(srcFile, outFile)


def CopyFilesToOutput(inputFolder, files, outputfolder = ""):
	for file in files:
		Show(file)
		CopyFileToOutput(inputFolder, file, outputfolder)	
