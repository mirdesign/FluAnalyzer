#!/bin/sh

#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#################################################################################################################################################################
# Use this script to:																																			#
# 1. Rename files based on an excel list																														#
#		TODO: description BLABKLABLABLA																															#
#                                                                                                                                                               #
# Expected Filename strategy in output: [__CUSTOM__]__f#[__CUSTOM__]-[segment]_[AA/NA].ext" e.g.                                                                #
#   [__CUSTOM__] = Custom text, ignored by the pipeline                                                                                                         #    
#   [segment] = Flu segment (i.e. PB2, HA, etc.)                                                                                                                #
# Keep the "__f#" (e.g. __f3) in the filename for automatic detection of the reading frame in subsequent steps of the pipeline                                  #
# Keep the -[segment]_[AA/NA].ext (e.g. -HA_AA.txt or -MP_NA.txt) unaltered for detection of segment and type in subsequent steps                               #
#                                                                                                                                                               #
#################################################################################################################################################################


#### INVARIABLE PART OF THE SCRIPT -- DO NOT CHANGE ANYTHING BELOW ####
. ./_system/shell/global.sh

#### Load config setting - Please change config to your needs ####
. ./_system/shell/config_main.sh

echo "###########################################"
echo "#                                         #"
echo "# Step 1. Housekeeping in file naming.    #"
echo "#                                         #"
echo "###########################################"
echo "Settings:"
echo -e "\ttext in filename to replace: ${string_to_replace}"
echo -e "\tnew text: ${new_string}"
echo -e "\tinput directory: ${DIR_working}"
echo

echo "Create backup of your original files:"
echo -e "\t${DIR_backup}/${BackupFile_Original_Sourcefiles}"
mkdir ${DIR_backup}
cd ${DIR_working}
zip -qr $BackupFile_Original_Sourcefiles *.* -x *.zip
mv $BackupFile_Original_Sourcefiles $DIR_backup

# Some housekeeping since the namings of the output files of the previous tool didn't comply with our pipeline 
echo
echo "Perform housekeeping of file names:"
echo -e "\t-Add an extra underscore in front of '_f' "
rename 's/([^\_]+)(\_f)/$1_$2/g' *.txt -vn

echo -e "\t-Replace T for URT"
rename 's/\_T_/_URT_/g' *.txt -vn

echo -e "\t-Replace N for URT"
rename 's/\_N_/_URT_/g' *.txt -vn

echo -e "\t-Replace parts of filename"
echo

python ${DIR_scripts_python}/replace_part_of_filename.py ${DIR_working} ${string_to_replace} ${new_string}
if [ $? -ne 0 ]; then >&2 echo "Python script failed. We quit."; exit; fi
echo
echo "Housekeeping ready!"
echo

. $DIR_scripts_shell/Auto_Proceed.sh $ProceedAfterShellEndOnSucces $shell_script_step_2
