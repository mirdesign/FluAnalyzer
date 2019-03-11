#!/bin/sh

#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################

#### VARIABLES CONFIG FILE AT PATH BELOW -- UPDATE BEFORE EACH ANALYSIS ####
# Syntax: ". [folder name]" (NOTE: Start with ".")

. ./_system/shell/global.sh
. ./_system/shell/config_main.sh


#### INVARIABLE PART OF THE SCRIPT -- DO NOT CHANGE ANYTHING BELOW ####
_input_directory_vd_gpi=''$DIR_working'/input_vd_gpi'
_folder_output_vd_gpi_mv_subfolder_name='output_step3'
_folder_output_vd_gpi_mv=''$DIR_working'/'$_folder_output_vd_gpi_mv_subfolder_name''
_folder_experiments=''$_folder_output_vd_gpi_mv'/experiments'
_experiment_timestamp_name=''$experiment_name'_'$(date +%s)''
_file_experiment='experiment_'${experiment_name}'.txt'
_filepath_experiment=''${_input_directory_vd_gpi}'/'${_file_experiment}''
_folder_experiment=''$_folder_experiments'/'$_experiment_timestamp_name''
_folder_variant_detector_output_mv_AA_folderName='out_mv_AA'
_folder_variant_detector_output_mv_nt_folderName='out_mv_nt'
_folder_variant_detector_output_mv_AA=''$_folder_experiment'/'$_folder_variant_detector_output_mv_AA_folderName''
_folder_variant_detector_output_mv_nt=''$_folder_experiment'/'$_folder_variant_detector_output_mv_nt_folderName''
_folder_variant_detector_output_mv_AA_primer_incl=''$_folder_variant_detector_output_mv_AA'/pr_inc'
_folder_variant_detector_output_mv_AA_primer_excl=''$_folder_variant_detector_output_mv_AA'/pr_ex'
_folder_variant_detector_output_mv_nt_primer_incl=''$_folder_variant_detector_output_mv_nt'/pr_inc'
_folder_variant_detector_output_mv_nt_primer_excl=''$_folder_variant_detector_output_mv_nt'/pr_ex'

_folder_variant_detector_output_mv_AA_primer_incl_overview=''$_folder_variant_detector_output_mv_AA_primer_incl'/overview'
_folder_variant_detector_output_mv_AA_primer_excl_overview=''$_folder_variant_detector_output_mv_AA_primer_excl'/overview'
_folder_variant_detector_output_mv_nt_primer_incl_overview=''$_folder_variant_detector_output_mv_nt_primer_incl'/overview'
_folder_variant_detector_output_mv_nt_primer_excl_overview=''$_folder_variant_detector_output_mv_nt_primer_excl'/overview'

_folder_variant_detector_output_AA_folderName='out_vd_AA'
_folder_variant_detector_output_nt_folderName='out_vd_nt'
_folder_variant_detector_output_AA=''$_folder_experiment'/'$_folder_variant_detector_output_AA_folderName''
_folder_variant_detector_output_nt=''$_folder_experiment'/'$_folder_variant_detector_output_nt_folderName''
_folder_variant_detector_input_AA=''${_folder_experiment}'/input_variant_detector_AA'
_folder_variant_detector_input_nt=''${_folder_experiment}'/input_variant_detector_nt'
_folder_output_positions_of_interest=''$_folder_experiment'/positions_of_interest'
_folder_output_r_visuals=''$_folder_experiment'/visuals'
_folder_variant_detector_output_mv_visuals=''$_folder_output_r_visuals'/mv'
_file_variant_detector_output_mv_visual_AA_incl='.AA.incl.pdf'
_file_variant_detector_output_mv_visual_AA_excl='.AA.excl.pdf'
_file_variant_detector_output_mv_visual_nt_incl='.nt.incl.pdf'
_file_variant_detector_output_mv_visual_nt_excl='.nt.excl.pdf'


####
##
#	Introduction
##
####
echo
echo
echo "#########################################"
echo -e "#                                     \t#"
echo -e "#\t Variant detector               #"
echo -e "#\t Minority variant detection     #"
echo -e "#\t Gather positions of interest   #"
echo -e "#\t Visualisation                  #"
echo -e "#                                     \t#"
echo -e "# Version: 1.0                        \t#"
echo "#########################################"
echo
cd ${_input_directory_vd_gpi}
if [ $(ls -l | grep -c $_file_experiment) -eq 0 ]; then
  echo
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "The experiment file does not exist in the current experiment directory!"
  echo "Terminating shell script to protect any previous output."
  echo "Current experiment directory: '$_input_directory_vd_gpi'"
  echo "Search for experiment file: '$_file_experiment'"
  echo "Did you remove 'experiment_' and the extention (i.e. '.txt') in the configuration setting?"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo
  exit 1
fi
echo "Used settings:"
if [ "${run_script_in_parallel}" != 1 ]; then
	echo -e "\tRun in parallel: No"
else
	echo -e "\tRun in parallel: Yes"
fi
echo -e "\tBasic analysis files: $DIR_basic_analysis_files"
echo -e "\tExperiment name: $experiment_name"
echo -e "\tExperiment output folder: $_folder_experiment"
echo -e "\tPrimer list: $primer_list_AA"
echo -e "\tPython scripts: $DIR_scripts_python"
echo -e "\tSamples to analyse:"  
cat $_filepath_experiment | while read line; do
    echo -e "\t\t- $line"
done
echo

####
##
#	Step 3.1 Variant Detector
##
####

# Add a blank line to the end of $_file_experiment. The last iteration of the file is skipped if no blank line at the end of the file 
sed -i '$ a \' $_filepath_experiment

echo
echo "#####################################################"
echo
echo "STEP 3.1: Variant detection"
echo
echo "Create directories..."

cd ${DIR_working}
if [ $(ls -l | grep -c ${_folder_output_vd_gpi_mv_subfolder_name}) -eq 0 ]; then
	mkdir ${_folder_output_vd_gpi_mv}
fi

cd ${_folder_output_vd_gpi_mv}
if [ $(ls -l | grep -c "experiments") -eq 0 ]; then
	mkdir ${_folder_experiments}
fi

echo "Validate experiment location..."
cd ${_folder_experiments}
if [ $(ls -l | grep -c ${_experiment_timestamp_name}) -gt 0 ]; then
  echo
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "The experiment folder exist in the current experiment directory,"
  echo "therefore we assume a process has been run before."
  echo "Terminating shell script to protect any previous output."
  echo "Current experiment directory: '$_folder_experiment'"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo
  exit 1
fi

mkdir ${_folder_experiment}
mkdir ${_folder_variant_detector_output_mv_AA}
mkdir ${_folder_variant_detector_output_mv_nt}
mkdir ${_folder_variant_detector_output_mv_AA_primer_incl}
mkdir ${_folder_variant_detector_output_mv_AA_primer_excl}
mkdir ${_folder_variant_detector_output_mv_nt_primer_incl}
mkdir ${_folder_variant_detector_output_mv_nt_primer_excl}
mkdir ${_folder_variant_detector_input_AA}
mkdir ${_folder_variant_detector_input_nt}
mkdir ${_folder_variant_detector_output_AA}
mkdir ${_folder_variant_detector_output_nt}
mkdir ${_folder_output_r_visuals}
mkdir ${_folder_variant_detector_output_mv_visuals}
mkdir ${_folder_variant_detector_output_mv_AA_primer_incl_overview}
mkdir ${_folder_variant_detector_output_mv_AA_primer_excl_overview}
mkdir ${_folder_variant_detector_output_mv_nt_primer_incl_overview}
mkdir ${_folder_variant_detector_output_mv_nt_primer_excl_overview}



echo "Copy data for current experiment..."
while read file; do 
cd ${DIR_basic_analysis_files}/AA-files/${file}/$DIRNAME_Output_Step2_Protein
if [ $(ls -l | grep -c "txt") -gt 0 ]; then	
	cp *.txt ${_folder_variant_detector_input_AA} &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
fi
cd ${DIR_basic_analysis_files}/nt-files/${file}/$DIRNAME_Output_Step2_Gene
if [ $(ls -l | grep -c "txt") -gt 0 ]; then	
	cp *.txt ${_folder_variant_detector_input_nt} &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
fi
done < $_filepath_experiment


# Variant Detector AA
if [ "${run_script_in_parallel}" = 1 ]; then
echo
echo -e "Execute ${Gre}VARIANT DETECTOR${RCol} for AA-files and nt-files in parallel"
echo -e $Gre
python ${DIR_scripts_python}/variant_detector.py ${_folder_variant_detector_input_AA} ${_folder_variant_detector_output_AA} ${primer_list_AA} ${experiment_name} ${Minimal_Value_RowCoverage} ${VarCovOK_Minimal_Value_Count} ${VarCovOK_Minimal_Value_Percentage} ${Minimal_Value_Percentage_Change} ${Minimal_Value_Row_Var_Perc} ${Minimal_Value_Percentage_Difference} ${Threshold_Value_Percentage_Difference} ${Threshold_Value_Max_Reached} ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time}  ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location} &
else
echo
echo -e "Execute ${Gre}VARIANT DETECTOR${RCol} for AA-files"
echo -e $Gre
python ${DIR_scripts_python}/variant_detector.py ${_folder_variant_detector_input_AA} ${_folder_variant_detector_output_AA} ${primer_list_AA} ${experiment_name} ${Minimal_Value_RowCoverage} ${VarCovOK_Minimal_Value_Count} ${VarCovOK_Minimal_Value_Percentage} ${Minimal_Value_Percentage_Change} ${Minimal_Value_Row_Var_Perc} ${Minimal_Value_Percentage_Difference} ${Threshold_Value_Percentage_Difference} ${Threshold_Value_Max_Reached} ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time}  ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location}
echo -e $RCol
if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi
fi

# Variant Detector nt
if [ "${run_script_in_parallel}" = 1 ]; then
python ${DIR_scripts_python}/variant_detector.py ${_folder_variant_detector_input_nt} ${_folder_variant_detector_output_nt} ${primer_list_nt} ${experiment_name} ${Minimal_Value_RowCoverage} ${VarCovOK_Minimal_Value_Count} ${VarCovOK_Minimal_Value_Percentage} ${Minimal_Value_Percentage_Change} ${Minimal_Value_Row_Var_Perc} ${Minimal_Value_Percentage_Difference} ${Threshold_Value_Percentage_Difference} ${Threshold_Value_Max_Reached} ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time}  ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location} &
else
echo
echo -e "Execute ${Blu}VARIANT DETECTOR${RCol} for nt-files"
echo -e $Blu
python ${DIR_scripts_python}/variant_detector.py ${_folder_variant_detector_input_nt} ${_folder_variant_detector_output_nt} ${primer_list_nt} ${experiment_name} ${Minimal_Value_RowCoverage} ${VarCovOK_Minimal_Value_Count} ${VarCovOK_Minimal_Value_Percentage} ${Minimal_Value_Percentage_Change} ${Minimal_Value_Row_Var_Perc} ${Minimal_Value_Percentage_Difference} ${Threshold_Value_Percentage_Difference} ${Threshold_Value_Max_Reached} ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Time}  ${RowMultiVarCovOK_Minimal_Value_Amount_TimePoints_VarRowCoverageOk_Location}
if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi &
fi
echo -e $RCol
wait

# Archive
zip -qrm ${_folder_experiment}/input_variant_detector.zip $_folder_variant_detector_input_AA
zip -qrm ${_folder_experiment}/input_variant_detector.zip $_folder_variant_detector_input_nt


####
##
#	Step 3.2 Minority variant
##
####

echo
echo "################################################################"
echo "STEP 3.2: Identify and report AA and nt positions of interest"
echo 
echo

## Minority variant detection
echo -e "Execute ${Yel}MINORITY VARIANT CALCULATOR${RCol} for AminoAcids and Nucleotides, incl/excl primer regions"
echo
echo -e "\t${RCol}* Including primer regions for AminoAcids${Yel}"
python ${DIR_scripts_python}/minority_variant_calculator.py ${_folder_variant_detector_output_AA}/primer_incl ${_folder_variant_detector_output_mv_AA_primer_incl} ${_folder_variant_detector_output_mv_AA_primer_incl_overview} &
if [ "${run_script_in_parallel}" != 1 ]; then wait; if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi; fi

echo -e "\t${RCol}* Excluding primer regions for AminoAcids${Yel}"
python ${DIR_scripts_python}/minority_variant_calculator.py ${_folder_variant_detector_output_AA}/primer_excl ${_folder_variant_detector_output_mv_AA_primer_excl} ${_folder_variant_detector_output_mv_AA_primer_excl_overview} &
if [ "${run_script_in_parallel}" != 1 ]; then wait; if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi; fi

echo -e "\t${RCol}* Including primer regions for Nucleotides${Yel}"
python ${DIR_scripts_python}/minority_variant_calculator.py ${_folder_variant_detector_output_nt}/primer_incl ${_folder_variant_detector_output_mv_nt_primer_incl} ${_folder_variant_detector_output_mv_nt_primer_incl_overview} &
if [ "${run_script_in_parallel}" != 1 ]; then wait; if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi; fi

echo -e "\t${RCol}* Excluding primer regions for Nucleotides${Yel}"
python ${DIR_scripts_python}/minority_variant_calculator.py ${_folder_variant_detector_output_nt}/primer_excl ${_folder_variant_detector_output_mv_nt_primer_excl} ${_folder_variant_detector_output_mv_nt_primer_excl_overview} &
wait
if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi
echo -e $RCol


## Minority variant visualisation
echo
echo "Build Minority Variant visuals for AminoAcids and Nucleotides, incl/excl primer regions"
echo
echo -e "\t${RCol}* Including primer regions for AminoAcids${Yel}"
Rscript ${DIR_scripts_python}/minority_variant_visual.R ${_folder_variant_detector_output_mv_AA_primer_incl} ${_folder_variant_detector_output_mv_visuals} ${_file_variant_detector_output_mv_visual_AA_incl} &
if [ "${run_script_in_parallel}" != 1 ]; then wait; fi

echo -e "\t${RCol}* Excluding primer regions for AminoAcids${Yel}"
Rscript ${DIR_scripts_python}/minority_variant_visual.R ${_folder_variant_detector_output_mv_AA_primer_excl} ${_folder_variant_detector_output_mv_visuals} ${_file_variant_detector_output_mv_visual_AA_excl} &
if [ "${run_script_in_parallel}" != 1 ]; then wait; fi

echo -e "\t${RCol}* Including primer regions for Nucleotides${Yel}"
Rscript ${DIR_scripts_python}/minority_variant_visual.R ${_folder_variant_detector_output_mv_nt_primer_incl} ${_folder_variant_detector_output_mv_visuals} ${_file_variant_detector_output_mv_visual_nt_incl} &
if [ "${run_script_in_parallel}" != 1 ]; then wait; fi

echo -e "\t${RCol}* Excluding primer regions for Nucleotides${Yel}"
Rscript ${DIR_scripts_python}/minority_variant_visual.R ${_folder_variant_detector_output_mv_nt_primer_excl} ${_folder_variant_detector_output_mv_visuals} ${_file_variant_detector_output_mv_visual_nt_excl} &
wait
echo -e $RCol


####
##
#	Step 3.3 Gather positions of interest
##
####
echo
echo "##############################################################"
echo "STEP 3.3: Identify and report AA and nt positions of interest"
echo 

echo "Create output directory..."
mkdir $_folder_output_positions_of_interest

# Provide primer_incl folder as input to generate two output files one with and one without primer positions listed
echo
echo -e "Execute ${Cya}GATHER POSITIONS OF INTEREST${RCol} for AA-files and nt-files"
echo -e $Cya
python ${DIR_scripts_python}/gather_positions_of_interest.py ${_folder_variant_detector_output_nt}/primer_incl ${_folder_variant_detector_output_AA}/primer_incl ${_folder_output_positions_of_interest} ${experiment_name}
if [ $? -ne 0 ]; then >&2 echo -e "${BRed}Python script failed. We quit."; exit; fi &
echo -e $RCol

echo
echo "##############################################################"
echo "STEP 3.4: Create visuals"
echo 
echo

# extract information out of the sample name
_compareGroupName=''$(echo ${experiment_name} | cut -d'_' -f 3)'' # H3N2_L1284_T91_NK.TA -> T91
_compareGroupName_FirstChar=''$(echo ${_compareGroupName} | cut -c1-1)'' # T91 -> T
_compareGroupName_RemainingChars=''$(echo ${_compareGroupName} | cut -c2-)'' # T91 -> 91
_compareGroup_IsTimePoint="False"

if [ "${_compareGroupName_FirstChar}" = "T" ]; then
	if echo $_compareGroupName_RemainingChars | egrep -q '^[0-9]+$'; then
		_compareGroup_IsTimePoint="True"
	fi
fi

if [ "${_compareGroup_IsTimePoint}" = "True" ]; then
	echo "Execute R script (Sample based) to create visuals for:"
	echo -e "\t* AminoAcids"
	echo -e $Pur
	Rscript ${DIR_scripts_python}/variant_detector_visual_locations.R ${_folder_output_positions_of_interest}/${experiment_name}.AA.csv ${_folder_output_r_visuals}/${experiment_name}.AA.pdf &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
	
	echo -e "\t${RCol}* AminoAcids - Non primer positions${Pur}"
	Rscript ${DIR_scripts_python}/variant_detector_visual_locations.R ${_folder_output_positions_of_interest}/${experiment_name}.AA.NonPrimersOnly.csv ${_folder_output_r_visuals}/${experiment_name}.AA.NonPrimersOnly.pdf &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
	
	echo -e "\t${RCol}* Nucleotides${Pur}"
	Rscript ${DIR_scripts_python}/variant_detector_visual_locations.R ${_folder_output_positions_of_interest}/${experiment_name}.NT.csv ${_folder_output_r_visuals}/${experiment_name}.NT.pdf &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
	
	echo -e "\t${RCol}* Nucleotides - Non primer positions${Pur}" 
	Rscript ${DIR_scripts_python}/variant_detector_visual_locations.R ${_folder_output_positions_of_interest}/${experiment_name}.NT.NonPrimersOnly.csv ${_folder_output_r_visuals}/${experiment_name}.NT.NonPrimersOnly.pdf &
	wait
else
	echo "Execute R script (Timepoints based) to create visuals for:"
	echo -e "\t* AminoAcids${Pur}"
	Rscript ${DIR_scripts_python}/variant_detector_visual_timepoints.R ${_folder_output_positions_of_interest}/${experiment_name}.AA.csv ${_folder_output_r_visuals}/${experiment_name}.AA.pdf &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
	
	echo -e "\t${RCol}* AminoAcids - Non primer positions${Pur}"
	Rscript ${DIR_scripts_python}/variant_detector_visual_timepoints.R ${_folder_output_positions_of_interest}/${experiment_name}.AA.NonPrimersOnly.csv ${_folder_output_r_visuals}/${experiment_name}.AA.NonPrimersOnly.pdf &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
	
	echo -e "\t${RCol}* Nucleotides${Pur}"
	Rscript ${DIR_scripts_python}/variant_detector_visual_timepoints.R ${_folder_output_positions_of_interest}/${experiment_name}.NT.csv ${_folder_output_r_visuals}/${experiment_name}.NT.pdf &
	if [ "${run_script_in_parallel}" != 1 ]; then wait; fi
	
	echo -e "\t${RCol}* Nucleotides - Non primer positions${Pur}" 
	Rscript ${DIR_scripts_python}/variant_detector_visual_timepoints.R ${_folder_output_positions_of_interest}/${experiment_name}.NT.NonPrimersOnly.csv ${_folder_output_r_visuals}/${experiment_name}.NT.NonPrimersOnly.pdf &
	wait
fi
wait
echo -e "\t${RCol}* Minority Variants - Nucleotides${Blu}" 
Rscript ${DIR_scripts_python}/minority_variant_visual_overview.R ${_folder_variant_detector_output_mv_nt_primer_incl_overview}/overview.csv ${_folder_variant_detector_output_mv_visuals}/${experiment_name}.nt.mv.pdf &
if [ "${run_script_in_parallel}" != 1 ]; then
wait
fi

echo -e "\t${RCol}* Minority Variants - Nucleotides - Non primer positions${Blu}" 
Rscript ${DIR_scripts_python}/minority_variant_visual_overview.R ${_folder_variant_detector_output_mv_nt_primer_excl_overview}/overview.csv ${_folder_variant_detector_output_mv_visuals}/${experiment_name}.nt.mv.NonPrimersOnly.pdf &
if [ "${run_script_in_parallel}" != 1 ]; then
wait
fi

echo -e "\t${RCol}* Minority Variants - AminoAcids${Blu}" 
Rscript ${DIR_scripts_python}/minority_variant_visual_overview.R ${_folder_variant_detector_output_mv_AA_primer_incl_overview}/overview.csv ${_folder_variant_detector_output_mv_visuals}/${experiment_name}.AA.mv.pdf &
if [ "${run_script_in_parallel}" != 1 ]; then
wait
fi

echo -e "\t${RCol}* Minority Variants - AminoAcids - Non primer positions${Blu}" 
Rscript ${DIR_scripts_python}/minority_variant_visual_overview.R ${_folder_variant_detector_output_mv_AA_primer_excl_overview}/overview.csv ${_folder_variant_detector_output_mv_visuals}/${experiment_name}.AA.mv.NonPrimersOnly.pdf &
wait
echo -e $RCol

echo "Archive working data..."
zip -qrm $_folder_experiment/output_vd_mv.zip ${_folder_variant_detector_output_mv_AA}
zip -qrm $_folder_experiment/output_vd_mv.zip ${_folder_variant_detector_output_mv_nt}
zip -qrm $_folder_experiment/output_vd_mv.zip ${_folder_variant_detector_output_AA}
zip -qrm $_folder_experiment/output_vd_mv.zip ${_folder_variant_detector_output_nt}
wait

echo "Done"