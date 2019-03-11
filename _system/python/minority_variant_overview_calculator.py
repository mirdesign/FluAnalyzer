#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# Minority Variant Overview Calculator                      #
#############################################################
import locale
import sys
import collections
import os
import glob
import re
import math
from shutil import copyfile

# Load common shared variables and methods
execfile(os.path.join(os.path.dirname(__file__),"_global_constants.py"))

inputFolders = []
outputFolder = ""

# Constants
__FILE_CONTENT_MINORITY_VARIANT = "PropMinVar (%)"
__FILE_CONTENT_MAJORITY_VARIANT = "PropMajVar (%)"
__FILE_CONTENT_HEADER_POSITION = "position"

__FILE_VARIANT_OVERVIEW_FILE_EXTENTION = "csv"
		
		
# Show some startup information
def ShowSplash():
	Show ("")
	Show ("--------------------------------")
	Show ("Tool: Minority Variant Overview Calculator")
	Show ("Take the output of Minority Variant Calculater and calculate the 'Minority variant frequency' over several experiments.")
	Show ("version: 1.0")
	Show ("Contact: m.pater@amc.nl | mpater@mirdesign.nl - Maarten Pater")
	Show ("--------------------------------")
	Show ("Usage: python minority_variant_overview_calculator.py folder_output {folders_with_output_minority_variant_calculator}")
	Show ("Please provide the output folders from of Minority Variant Calculator.")
	Show ("--------------------------------")
		
	
def mean(data):
   # "Return the sample arithmetic mean of data."
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/n # in Python 2 use sum(data)/float(n)

def _ss(data):
    #"Return sum of square deviations of sequence data."
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data):
    #"Calculates the population standard deviation."
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5
	
# Detect header locations
def GetColumnPositionFromHeader(line, SearchForHeader):
	pos = line.find(SearchForHeader)
			
	# Check if we found the headers
	if (pos == -1):
		QuitAndShowError([	("(GetColumnPositionFromHeader) " + SearchForHeader + " value could not be detected in the headerline."),
							("Provided: " + line),
							("Searched for: " + SearchForHeader)])
							
	return line[:pos].count(__INPUT_ROW_SEPERATOR)


# Build CSV header
def buildCSVHeader(file):
	file.write("position")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("mean")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("SEM")
	file.write(__OUTPUT_ROW_SEPERATOR)
	file.write("\n")
	
# Save CSV data to file	
# minorityVariantsPerSegment: segment, countHeader, pos, frequency
def saveCSV(MEANs, SEMs, numPositions):
	for fileName in MEANs:
		file = open(outputFolder + "/" + fileName, "w")
		file.write("sep=" + __OUTPUT_ROW_SEPERATOR + "\n")
		buildCSVHeader(file)
		
		for pos in range(1, numPositions[fileName] + 1):
		
			file.write(str(pos))
			file.write(__OUTPUT_ROW_SEPERATOR)

			if(pos in MEANs[fileName]):
				file.write(str(MEANs[fileName][pos]))
				file.write(__OUTPUT_ROW_SEPERATOR)
				file.write(str(SEMs[fileName][pos]))
			else:
				file.write("-")
				file.write(__OUTPUT_ROW_SEPERATOR)
				file.write("-")
					
			file.write("\n")
					
		file.close()



def ReadDataFromFile(file, fileFolder):
	
	skippedFirstLine = False
	headerLine = True
	
	colProteinPosition = -1
	colMinVar = -1
	numPositions = -1
	
	# pos, minVar
	retVal = {}
	
	for line in open(fileFolder + "/" + file):
		if(not skippedFirstLine): # Skip header
			skippedFirstLine = True
			continue
		
		if(headerLine):		
			
			colProteinPosition = GetColumnPositionFromHeader(line, __FILE_CONTENT_HEADER_POSITION)
			colMinVar = GetColumnPositionFromHeader(line, __FILE_CONTENT_MINORITY_VARIANT)
			
			
			headerLine = False
			continue
		
	
		cProteinPosition = int(GetDataFromFileContent(line, colProteinPosition))
		minVariant = GetDataFromFileContent(line, colMinVar)
		
		numPositions = int(cProteinPosition)
		
		if(not(minVariant == "-")):
			retVal[cProteinPosition] = float(minVariant)
			
	return retVal, numPositions
			
			
# rawDataFromFile format: # folder, pos, minVariant
# return: countHeader, pos, minorityVariant 
def GenerateMinorityVariantOverview (rawDataFromFile):

	means = {}
	SEM = {}
	positionValues = {}
	# Translate multiple folders into one position
	for folder in rawDataFromFile:
		for pos in rawDataFromFile[folder]:
			if not(pos in positionValues):
				positionValues[pos] = []
			positionValues[pos].append(rawDataFromFile[folder][pos])
	
	# Calculate means
	for pos in positionValues:
		numValues = len(positionValues[pos])
		if(numValues > 1):
			means[pos] = sum(positionValues[pos]) / numValues
			SEM[pos] = pstdev(positionValues[pos]) / math.sqrt(numValues)

	return means, SEM
		
# Process found files
# return: segment, countHeader, pos, frequency
def GetOverviewPerSegment(Files, fileFolders):
	if(len(Files) == 0):
		Show ("No files to process. Quit program")
		sys.exit(1)
			
	means = {}
	SEMs = {}
	numPositions = {}
		
	for file in Files:
		Show("Analysing files: %s" % file)
		rawDataFromFile = {}
		
		filesToCompare = 0
		
		for folder in fileFolders:
			if(os.path.isfile(folder + "/" + file)):
				filesToCompare = filesToCompare + 1
		
		
		if(filesToCompare > 1):
			for folder in fileFolders:
				rawDataFromFile[folder], numPositions[file] = ReadDataFromFile(file, folder) # pos, minVariant
			
			means[file], SEMs[file] =  GenerateMinorityVariantOverview(rawDataFromFile) # pos, mean  | pos, SEM
				
	return means, SEMs, numPositions
			
# Write log file
def WriteLogFile(Files):	
	# todo: create proper log file
	file = open(outputFolder + "/" + IdentifiableExperimentName + ".log", "w")
	file.write("Log file for minority variant overview calculator\n")
	file.write("\n")
	file.write("SETTINGS\n")
	file.write("showMessages: " + ("TRUE" if showMessages else "FALSE")+ "\n")
	file.write("deepDebug: " + ("TRUE" if deepDebug else "FALSE") + "\n")
	file.write("tabularDivide: " + tabularDivide + "\n") 
			
	file.write("inputFolder: " + inputFolder + "\n")
	file.write("outputFolder: " + outputFolder + "\n")
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
			
			
	outputFolder = sys.argv[1]
	for i in range(2, args):
		inputFolders.append(sys.argv[i])
	
		
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
	Show("Searching for minority variant files in the first folder...")
	Show("Folder: " + inputFolders[0])
	Files = GetFileDirectPaths(inputFolders[0])
	Files = RemoveNonSegmentedOverviewFiles(Files)


	Show("")
	Show ("Create minority variant overview...")
	MEANs, SEMs, numPositions = GetOverviewPerSegment(Files, inputFolders) # file, pos, mean | file, pos, SEM
		
	# (key, value)
	# [segment, variants]
	Show("Write output file...")
	saveCSV(MEANs, SEMs, numPositions)

	Show("Ready!") 
		
Main()