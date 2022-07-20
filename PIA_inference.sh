#!/bin/bash

#--------------------
# Author: Hesham ElAbd 
# Brief: A wrapper for the PIA-inference engine used to run the backend part of the server 
# Version: 0.3
# Copyright: Instuiute of Clinical Molecular Biology (UKSH), Kiel, Germany 2022. 
#----------------------

# Define a usage parameters 
PROGRAM_NAME=$0
VERSION='0.2'

# define A function with usage information 
function usage()
{
    echo "$PROGRAM_NAME Version:$VERSION:"
    echo "usage: $PROGRAM_NAME -d dir_name  -s standard_input -g genotype_table -f fasta_file -v vcf_file -a list_allele -t tissue_name -m model_index -w window_size -z step_size"
    echo 
    echo "    -d : The dir_name which is the path to the input directory where all files are located"    
    echo "    -s : The standard_input which is a table made of two or three columns, this table is fed directory to the inference engine."
    echo "    -g : The genotype_table, which is a tablge containing genotype information for generating personalized protein sequence from input sequnece,\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -f : The fasta_file, which is a fasta file containing protein sequences, the file shall contain the same protein ids as defined in the genotype-table.\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -v:  The vcf_file: the input to a *PHASED* and consequently called file that will be used for generating input protein sequences from each patient, 
        for more information check VCF2prot. This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -a: The list_allele: a file containing a list of alleles for running predictions on the pipeline."
    echo "    -t: tissue_name, the name of the tissue, only used for making predictions using PIA-M models."
    echo "    -m: The model_index: The index of the model, can be only, 1, 2, 3 or 4, othervalues are not allowed."
    echo "    -w: The window_size: used for generating peptides using a sliding window approach, define the windows size.\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -z: The step_size: used for generating peptides using a sliding window approach, define the step size.\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    
    echo
}
## Declare the variable into some default values 
#-----------------------------------------------
dir_name=-1; standard_input=-1; genotype_table=-1; 
fasta_file=-1; vcf_file=-1; list_allele=-1; 
model_index=-1; window_size=-1; step_size=-1;
tissue_name=-1; 

## If the program has been called without any parameter we print the help message and then quit the execution  
#------------------------------------------------------------------------------------------------------------
if [ "$#" -eq 1 ]; 
then 
    usage # print the usage help message 
    exit 0
fi 

## Parse command line arguments
#------------------------------
while getopts "d:s:g:f:v:a:m:w:z:t:h" flag;
do 
    case ${flag} in 
        h) usage
        exit 0;;
        d) dir_name=${OPTARG};; 
        s) standard_input=${OPTARG};;
        g) genotype_table=${OPTARG};;
        f) fasta_file=${OPTARG};;
        v) vcf_file=${OPTARG};;
        a) list_allele=${OPTARG};; 
        m) model_index=${OPTARG};; 
        w) window_size=${OPTARG};;
        t) tissue_name=${OPTARG};;
        z) step_size=${OPTARG};;
        \?) echo "I can not understand $OPTARG, it not one of the supported parameters. $PROGRAM_NAME can be called as follow " >& 2
            usage
            exit 1;;
        :) echo "Parameter $OPTARG. $PROGRAM_NAME can be called as follow " >& 2
            usage
            exit 1;;
    esac
done

## Parsing and checking the parameters  
#-------------------------------------


# 1. check the input directory exists 
if [ $dir_name == -1 ]
then 
    echo "ERROR:: $(date -u): Input directory has not been provide, update the root directory and try again." 
    exit 1
elif [ ! -d "$dir_name" ] 
then
    echo "ERROR:: $(date -u): The provided input path root directory: $dir_name does not exists." 
    exit 1
fi

## Building the outut structure 
#------------------------------
ROOT_DIR=$dir_name

# Building the Prelude 
#--------------------------------------
mkdir "$ROOT_DIR/output"
LAST_ERROR_CODE=$0

# checking for error with creating the output directory
#------------------------------------------------------
if [[ $LAST_EXIT_CODE != 0 ]]
then 
    if [[ -d "$ROOT_DIR/output" ]]
    then
        rm -r $ROOT_DIR/output
        mkdir "$ROOT_DIR/output"
        LAST_ERROR_CODE=$0

        if [[ $LAST_EXIT_CODE != 0 ]]
        then 
            echo "ERROR: $(date -u): Failed with the following error: $LAST_EXIT_CODE while creating a directory to store the results" 
            exit $LAST_EXIT_CODE # return the same exit code as the failed process 
        fi 
    else
    echo "ERROR: $(date -u): Failed with the following error: $LAST_EXIT_CODE while creating a directory to store the results" 
    exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi 
fi 

mkdir "$ROOT_DIR/stat"
LAST_ERROR_CODE=$0

# hecking for error with creating the stat directory
#---------------------------------------------------
if [[ $LAST_EXIT_CODE != 0 ]]
then 
    if [[ -d "$ROOT_DIR/stat" ]]
    then
        rm -r "$ROOT_DIR/stat"
        mkdir "$ROOT_DIR/stat"
        LAST_ERROR_CODE=$0

        if [[ $LAST_EXIT_CODE != 0 ]]
        then 
            echo "ERROR: $(date -u): Failed with the following error: $LAST_EXIT_CODE while creating a directory to store the logs and progress lines" 
            exit $LAST_EXIT_CODE # return the same exit code as the failed process 
        fi 
    else
    echo "ERROR: $(date -u): Failed with the following error: $LAST_EXIT_CODE while creating a directory to store the logs and progress lines" 
    exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi  
fi 

# 2. check that the logic of input command line is valid and set th execution Route
#----------------------------------------------------------------------------------
if [ -f "$dir_name/$standard_input" ] # IF Standard input is defined then it has the highest precedence and we are executing using first execution route
then 
    ROUTE=1
elif [[ -f "$dir_name/$genotype_table" ]] && [[ -f "$dir_name/$fasta_file" ]] # if a fasta file and a genotypte table are defined then we are going down the second execution route 
then
    ROUTE=2
elif [[ -f "$dir_name/$vcf_file" ]] && [[ -f "$dir_name/$fasta_file" ]]
then 
    ROUTE=3
elif [[ -f "$dir_name/$fasta_file" ]]
then
    ROUTE=4
else
    echo "ERROR:: $(date -u): Unknow execution route, requirments for running any of the execution pipelines has not been provided, the requirment as follow" 2>> "$ROOT_DIR/output/run.log"
    echo "Execution pipeline 1: Using a standard input table --> provide the table using -s" 2>> "$ROOT_DIR/output/run.log"
    echo "Execution pipeline 2: Using a simple genotype table (-g) along with a Reference fasta file (-f) and a list of alleles (-a)" 2>> "$ROOT_DIR/output/run.log"
    echo "Execution pipeline 3: Using a VCF file (-v), fasta file (-f) along with a list of alleles (-a)"2>> "$ROOT_DIR/output/run.log"
    echo "Execution pipeline 4: Using a fasta file (-f) and a list of alleles (-a)" 2>> "$ROOT_DIR/output/run.log"
    echo "provide the input for a valid execution pipeline and try again." 2>> "$ROOT_DIR/output/run.log"
    exit -1
fi
# Log the execution route to the log sink
echo "INFO:: $(date -u): Execution is going to follow ROUTE $ROUTE" >> "$ROOT_DIR/output/run.log"

if [[ $ROUTE == 2 ]] || [[ $ROUTE == 3 ]] || [[ $ROUTE == 4 ]]
then 
    # Check the list of alleles
    if ! [[ -f "$dir_name/$list_allele" ]]
    then 
        echo "ERROR: $(date -u): Unmet requriment for running execution pipeline number $ROUTE. List of alleles has not been provide (-a), provide the list of alleles and try again." 2>> "$dir_name/run.log"
        exit -1;
    fi 
    # Check the window size and set the default value
    if [[ $model_index == -1 ]]
    then 
        echo "WARNING: $(date -u): The model index has not been selected ..., setting it a default value of 1." >> "$ROOT_DIR/stat/.warn"
        model_index=1
    fi 
    # Check the window size and set the default value
    if [[ $window_size == -1 ]]
    then 
        echo "WARNING: $(date -u): The window size has not been set ..., setting it a default value of 15." >> "$ROOT_DIR/stat/.warng"
        window_size=15
    fi 
    # Check the step size and set the default value
    if [[ $step_size == -1 ]]
    then
        echo "WARNING: $(date -u): The step size has not been set ..., setting it to a default value of 1." >> "$ROOT_DIR/stat/.warn"
        step_size=1
    fi
    # create a working directory for this route  
    mkdir -p "$dir_name/temp_work/creating_protein_sequences/" 
    if [ $? != 0 ]
    then 
        echo "Error Creating a directory for protein sequence generation at : $dir_name/temp_work/creating_protein_sequences/ . " 2>> "$ROOT_DIR/stat/run.error"
        exit 3
    fi 
fi 

# change the directory to the input file 
#---------------------------------------
cd $dir_name

# start executing the different pipelines
#-----------------------------------------
if [ $ROUTE == 1 ]
then 
    echo "INFO:: $(date -u): Starting the prediction with a file containing $( wc -l $standard_input)" >> "$ROOT_DIR/stat/run.error"
    prediction_engine_portal.py --base_dir "." \
    --input $standard_input --model_index $model_index\
    --output_path "./output/prediction_results.tsv"\
    --unmapped_results "./output/unmapped_results.tsv" >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error"
    echo "INFO:: $(date -u): Execution finished." >> "$ROOT_DIR/stat/run.error"
    rm -r "temp_work/creating_protein_sequences/" 
    exit 0

elif [ $ROUTE == 2 ]
then 
    # 1. Generating a protein sequence from mutation data
    Mut2Prot.py --input_table $genotype_table \
    --input_fasta $fasta_file --results_path "temp_work/creating_protein_sequences/fragmented_protein_sequences_precursor.tsv" \
    >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error" 

    let LAST_EXIT_CODE=$?
    if [ $LAST_EXIT_CODE != 0 ]
    then 
        echo "ERROR: $(date -u): Running Mut2Prot failed with the following error code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
        rm -r "temp_work/creating_protein_sequences/" 
        exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi 
    
    # 2. Remove the the 3rd column which contain the pheno
    cut -f1,2 temp_work/creating_protein_sequences/fragmented_protein_sequences_precursor.tsv > temp_work/creating_protein_sequences/fragmented_protein_sequences.tsv

    let LAST_EXIT_CODE=$?
    if [ $LAST_EXIT_CODE != 0 ]
    then 
        echo "ERROR: $(date -u): Cleaning the generated file failed with the following line of code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
        rm -r "temp_work/creating_protein_sequences/"
        exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi 

elif [ $ROUTE == 3 ]
then 
    # 1. make a directory to store the results
    if [ -d "temp_work/creating_protein_sequences/personalized_sequences" ]
    then 
        rm -r temp_work/creating_protein_sequences/personalized_sequences
        mkdir temp_work/creating_protein_sequences/personalized_sequences 
    else
        mkdir temp_work/creating_protein_sequences/personalized_sequences
    fi 

    # 2. Generate perfonalized protein sequences
        # 1.b running predictions 
    vcf2prot --vcf_file $vcf_file --fasta_ref $fasta_file -v -g mt -o "temp_work/creating_protein_sequences/personalized_sequences" \
        >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error" 
    
    let LAST_EXIT_CODE=$?
    if [ $LAST_EXIT_CODE != 0 ]
    then 
        echo "ERROR: $(date -u): Running vcf2prot failed with the following error code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
        rm -r "temp_work/creating_protein_sequences/"
        exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi      
    # 1. fragment the proteome of the input sample 
    fragmentor.py \
    --input_fasta_dir "temp_work/creating_protein_sequences/personalized_sequences" \
    --results_path "temp_work/creating_protein_sequences/fragmented_protein_sequences.tsv"\
    --window_size $window_size \
    --step_size $step_size \
    >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error"

    # 2. check that the fragmentation went with out problems  
    let LAST_EXIT_CODE=$?
    if [ $LAST_EXIT_CODE != 0 ]
    then 
        echo "ERROR: $(date -u): Running fragmentor with sample: $sample_proteome failed with the following error code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
        rm -r "temp_work/creating_protein_sequences/"
        exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi 

elif [ $ROUTE == 4 ]
then 
    # 1. Fragmenting the generated proteins 
    fragmentor.py --input_fasta_file $fasta_file \
    --results_path "temp_work/creating_protein_sequences/fragmented_protein_sequences.tsv" \
    --window_size $window_size --step_size $step_size \
    >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error"

    # Correct reading errors
    #-----------------------
    let LAST_EXIT_CODE=$?
    if [ $LAST_EXIT_CODE != 0 ]
    then 
        echo "ERROR: $(date -u): Running fragmentor.py failed with the following error code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
        rm -r "temp_work/creating_protein_sequences/"
        exit $LAST_EXIT_CODE # return the same exit code as the failed process 
    fi 
fi 
## Execute the common Route of each task
#----------------------------------------
# 3. Read the list of alleles and standrdize the names 
allele2standard.py --input_table $list_allele --output_table "temp_work/creating_protein_sequences/allele_names_stabdardized.tsv" \
    >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error"

# Correct reading errors
#-----------------------
let LAST_EXIT_CODE=$?
if [ $LAST_EXIT_CODE != 0 ]
then 
    echo "ERROR: $(date -u): Running Allele2standard failed with the following error code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
    rm -r "temp_work/creating_protein_sequences/"
    exit $LAST_EXIT_CODE # return the same exit code as the failed process 
fi 
#---------------------------------------------

# 4. Combine the alleles with the list of fragments
#--------------------------------------------------
Mixer.py --standardized_allele_name "temp_work/creating_protein_sequences/allele_names_stabdardized.tsv"\
    --input_peptide "temp_work/creating_protein_sequences/fragmented_protein_sequences.tsv" \
    --results_path "temp_work/creating_protein_sequences/input_to_the_prediction_engine.tsv" \
   >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error" 

# Checking for error and exit upon finding any
#---------------------------------------------
let LAST_EXIT_CODE=$?
if [ $LAST_EXIT_CODE != 0 ]
then 
    echo "ERROR: $(date -u): Running Mixer.py failed with the following error code: $LAST_EXIT_CODE" 2>> "$ROOT_DIR/stat/run.error"
    exit $LAST_EXIT_CODE # return the same exit code as the failed process 
fi 

# 5. run the predictions
#-----------------------
prediction_engine_portal.py --base_dir "." \
    --input "temp_work/creating_protein_sequences/input_to_the_prediction_engine.tsv"\
    --model_index $model_index \
    --output_path "../output/prediction_results.tsv" \
    --unmapped_results "../output/unmapped_results.tsv" \
    --tissue $tissue_name \
    >> "$ROOT_DIR/output/run.log" 2>> "$ROOT_DIR/stat/run.error"

# check for error and exit upon finding any
#------------------------------------------
let LAST_EXIT_CODE=$?
if [ $LAST_EXIT_CODE != 0 ]
then 
    echo "ERROR: $(date -u): Running prediction_engine_portal.py failed with the following error code: $LAST_EXIT_CODE " 2>> "$ROOT_DIR/stat/run.error"
    rm -r "temp_work/creating_protein_sequences/"
    exit $LAST_EXIT_CODE # return the same exit code as the failed process 
fi

# 7. clean up the temp directory
rm -r temp_work/creating_protein_sequences/
rm -r temp_work/
rm -r ./output/

echo "INFO:: $(date -u): Exectuion finish" >> "$ROOT_DIR/output/run.log"
# call PIA-inference wrapper engine 
#----------------------------------

