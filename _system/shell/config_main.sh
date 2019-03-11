#!/bin/sh
#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# This file concatenates basic setting values.              #
# Here you also load you specific project configuration     #
# file, which should be located in the '/config/' folder.   #
#############################################################

# dir where python and R scripts are located
DIR_scripts_python='/home/mpater/FluAnalyser/Pipeline/_system/python'
DIR_scripts_shell='/home/mpater/FluAnalyser/Pipeline/_system/shell'

#### INVARIABLE PART OF THE SCRIPT -- DO NOT CHANGE ANYTHING BELOW ####
. ./config/config_manual.sh

# Some global non-configurable variables
timestamp="$(date +%s)"

DIR_Start=$(pwd)
DIR_working=''${input_dir}'/'${current_analysis}''
DIR_log=''$DIR_working'/logfiles' # Path to log files directory
DIR_unique_sample_list=''$DIR_working'/__sample_list' # Path to DIR_unique_sample_list directory
FILE_Iterate_Sample_list_AA=''${DIR_unique_sample_list}'/filelist_reading_frames_AA-temp.isl'
FILE_Iterate_Sample_list_NT=''${DIR_unique_sample_list}'/filelist_reading_frames_nt-temp.isl'
DIR_input_vd_gpi=''$DIR_working'/input_vd_gpi'
DIR_basic_analysis_files=''${DIR_input_vd_gpi}'/basic_analysis_files'

DIR_output_rfs_pnc=''${DIR_working}'/output_rfs_pnc'
DIR_output_rfs_pnc_AA=''${DIR_output_rfs_pnc}'/AA-files'
DIR_output_rfs_pnc_NT=''${DIR_output_rfs_pnc}'/nt-files'
DIR_backup=''${DIR_working}'/backup'
FILE_Backup_UnusedFiles=''${DIR_backup}'/step2.non-selected_rf_files.zip'
FILE_Log_PNC_Warnings=''${DIR_log}'/step2.pnc_warnings.log'
FILE_Backup_RFS_PNC_Output=''${DIR_backup}'/step2.output_rfs_pnc.zip'
DIRNAME_Output_Step2_Protein="pnc_out_protein_files"
DIRNAME_Output_Step2_Gene="pnc_out_gene-segment_files"

primer_list_AA=''${primer_dir}'/'${primer_list_AA_file}''
primer_list_nt=''${primer_dir}'/'${primer_list_nt_file}''

BackupFile_Original_Sourcefiles="step1.original_files_backup__$timestamp.zip"

shell_script_step_1="1.file_housekeeping.sh"
shell_script_step_2="2.ReadingFrameSelector_and_ProteinNumberConverter.sh"

alias ls='ls' #we need ls to work normally