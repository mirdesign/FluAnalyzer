#!/bin/sh



### configuration file for sample 2015
# Name of the current analysis, which could be the name of the output directory of the previous step in the analysis (e.g. translate_sam_to_AA--job7F)
current_analysis='1205'

primer_list_AA_file='H3N2_primer_positions-AA.csv'
primer_list_nt_file='H3N2_primer_positions-nt.csv'

experiment_name="all"

# Path to directory holding the files to be renamed
input_dir='/home/mpater/FluAnalyser/data/patients'
primer_dir='/home/mpater/FluAnalyser/data/primer_position_identification_lists'



### configuration file for 1.NGS_file_renaming.sh
# part of the filename to replace 
string_to_replace='.2013-12-21-A.Stockholm.40.2013-EPI_ISL_155657-direct-A.H3N2-' 

## Run script in parallel mode?
## If yes, the output on screen will be messy but the running time will decrease significantly.
## 1 = yes
## 0 = no
run_script_in_parallel=1

# Do you want to auto continue the next shell script? 
#   "Continious" will automatically proceed.
#   "Ask" will ask before proceed.
#   "No" will not proceed.
#ProceedAfterShellEndOnSucces="No"
ProceedAfterShellEndOnSucces="Ask"
#ProceedAfterShellEndOnSucces="Continious"

# new text (Use "-" to remove a part of the filename as it replaces the full string ending in "-" with only "-") 
new_string="-"

#### INVARIABLE PART OF THE SCRIPT -- DO NOT CHANGE ANYTHING BELOW ####
# Check if no other config has been loaded already
# Please always copy this part when building a new configuration file
if [ -n "$__CONFIG_LOADED" ]; then
    echo "Please do not load multiple config files. Check config_main.sh for mistakes."
    exit
fi
__CONFIG_LOADED="${0}"