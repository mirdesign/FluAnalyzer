#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# Replace of remove part of a filename                      #
#############################################################
import os
import sys
print("")
print("------------------------------------------------------------------------")
print("Tool: replace / remove part of filename")
print("Version: 1.0")
print("Removes or replaces part of the filename of files in specified folder")
print("Usage: " + sys.argv[0] + ".py path_to_folder_with_files_to_rename string_to_be_replaced new_string")
print("To remove part of a filename use "" for new_string")
print("------------------------------------------------------------------------")
print("")

if len (sys.argv) != 4 :
    print("Incorrect number of arguments provided, please provide correct arguments.")
    sys.exit (1)

pathInput = sys.argv[1]		# dir with files
replaceFrom = sys.argv[2] 	# string to be replaced
replaceTo = sys.argv[3]		# new string. 

# Loop files
print ("Start renaming files.."),
debug=False
for root, subFolders, files in os.walk(pathInput):
	for file in files:
		oldName = os.path.join(root, file)
		newName = os.path.join(root, file.replace(replaceFrom, replaceTo))

		if(debug):
			print ("Consider: " + oldName)
		if(not oldName == newName):
			if(debug):
				print (" [Rename: " + oldName + " ==> " + newName + "]")
			os.rename(oldName, newName)
		elif(debug):
			print (" [Skip. Nothing to change]")

print("Ready!")