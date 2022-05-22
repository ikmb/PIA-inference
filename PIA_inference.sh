#!/bin/bash

# Author: Hesham ElAbd 
# Brief: A wrapper for the PIA-inference engine used to run the backend part of the server 
# Version: 0.2
# Copyright: Instuiute of Clinical Molecular Biology (UKSH), Kiel, Germany 2022. 

# Define a usage parameters 
PROGRAM_NAME=$0
VERSION='0.2'

# define A function with usage information 
function usage()
{
    echo "usage: $PROGRAM_NAME (Version:$VERSION) -d dir_name  -s standard_input -g genotype_table -f fasta_file -v vcf_file -a list_allele -t tissue_name -m model_index -w window_size -z step_size"
    echo 
    echo "    -d : The dir_name which is the path to the input directory where all files are located"    
    echo "    -s : The standard_input which is a table made of two or three columns, this table is feed directory to the inference engine."
    echo "    -g : The genotype_table, which is a tablge containing genotype information for generating personalized protein sequence from input sequnece,\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -f : The fasta_file, which is a fasta file containing protein sequences, the file shall contain the same protein ids as defined in the genotype-table.\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -v:  The vcf_file: the input to a *PHASED* and consequently called file that will be used for generating input protein sequences from each patient, 
        for more information check VCF2prot. This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
    echo "    -a: The list_allele: a file containing a list of alleles for running predictions on the pipeline.\
 This input can also be set to -1 which suggests to the pipeline that this step shall not be generated, i.e. do not generate personalized sequences."
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
## Parse command line arguments
#------------------
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

## parsing and checking the 
#------------------
# 1. check the input directory exists 
if [ $dir_name == -1 ]
then 
    echo "ERROR:: $(date -u): Input directory has not been provide, update the root directory and try again." 2>> "$dir_name/$( basename $dir_name).logs"
    exit 1
elif [ ! -d "$dir_name" ] 
then
    echo "ERROR:: $(date -u): The provided input path root directory: $dir_name does not exists." 2>> "$dir_name/$( basename $dir_name).logs"
    exit 1
fi

# 2. check that the logic of input command line is valid and set th execution Route
#----------------------------------------------------------------------------------
if [ -d $standard_input ] # IF Standard input is defined then it has the highest precedence and we are executing using first execution route
then 
    ROUTE=1
elif [[ -d $genotype_table ]] && [[ -d $fasta_file ]] # if a fasta file and a genotypte table are defined then we are going down the second execution route 
then
    ROUTE=2
elif [[ -d $vcf_file ]] && [[ -d $fasta_file ]]
then 
    ROUTE=3
elif [ -d $fasta_file ]
then
    ROUTE=4
else
    echo "ERROR:: $(date -u): Unknow execution route, requirments for running any of the execution pipelines has not been provided, the requirment as follow" 2>> "$dir_name/$( basename $dir_name).logs"
    echo "Execution pipeline 1: Using a standard input table --> provide the table using -s" 2>> "$dir_name/$( basename $dir_name).logs"
    echo "Execution pipeline 2: Using a simple genotype table (-g) along with a Reference fasta file (-f) and a list of alleles (-a)" 2>> "$dir_name/$( basename $dir_name).logs"
    echo "Execution pipeline 3: Using a VCF file (-v), fasta file (-f) along with a list of alleles (-a)"2>> "$dir_name/$( basename $dir_name).logs"
    echo "Execution pipeline 4: Using a fasta file (-f) and a list of alleles (-a)" 2>> "$dir_name/$( basename $dir_name).logs"
    echo "provide the input for a valid execution pipeline and try again." 2>> "$dir_name/$( basename $dir_name).logs"
    exit -1
fi
# Log the execution route to the log sink
echo "INFO:: $(date -u): Execution is going to follow ROUTE $ROUTE" >> "$dir_name/$( basename $dir_name).logs"

if [[ $ROUTE == 2 ]] || [[ $ROUTE == 3 ]] || [[ $ROUTE == 4 ]]
then 
    # Check the list of alleles
    if [ -d $list_allele]
    then 
        echo "ERROR: $(date -u): Unmet requriment for running execution pipeline number $ROUTE. List of alleles has not been provide (-a), provide the list of alleles and try again." 2>> "$dir_name/$( basename $dir_name).logs"
        exit -1;
    fi 
    # Check the window size and set the default value
    if [ $window_size == -1 ]
    then 
        echo "WARNING: $(date -u): The window size has not been set ..., setting it a default value of 15." >> "$dir_name/$( basename $dir_name).logs"
        window_size = 15
    fi 
    # Check the step size and set the default value
    if [ $step_size == -1 ]
    then
        echo "WARNING: $(date -u): The step size has not been set ..., setting it to a default value of 1." >> "$dir_name/$( basename $dir_name).logs"
        step_size = 1
    fi
    # create a working directory for this route  
    mkdir -p "$dir_name/temp_work/creating_protein_sequences/" 
    if [ echo $ != 0 ]
    then 
        echo "Error Creating a directory for protein sequence generation at : $dir_name/temp_work/creating_protein_sequences/ . " 2>> "$dir_name/$( basename $dir_name).err"
        exit 3
    fi 
fi 
# activate the conda environment 
#-------------------------------
source /home/helabd/miniconda3/bin/activate PIA_deployment

# export the path for calling PIA
#--------------------------------
export PATH=/home/helabd/PIA:$PATH

# change the directory to the input file 
#---------------------------------------
cd dir_name

# start executing the different pipelines
#-----------------------------------------
if [ $ROUTE == 1 ]
then 
    echo "INFO:: $(date -u): Starting the prediction with a file containing $( wc -l $standard_input)" >> "$dir_name/$( basename $dir_name).err"
    prediction_engine_portal.py --base_dir $dir_name \
    --input $standard_input --model_index $model_index\
    --output_path "output/prediction_results.tsv"\
    --unmapped_results "output/unmapped_results.tsv" >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err"
    echo "INFO:: $(date -u): Execution finished." >> "$dir_name/$( basename $dir_name).err"
    exit 0
elif [ $ROUTE == 2 ]
then 
    # 1. Generating a protein sequence from mutation data
    Mut2Prot.py --genetic_table $genotype_table \
    --fasta_file $fasta_file --output_dir "$dir_name/temp_work/creating_protein_sequences/personalized_protein_sequences.fasta" \
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err" 
        
    # 2. Fragment input collection of Fragments 
    fragmentor.py --input_fasta_file "$dir_name/temp_work/creating_protein_sequences/personalized_protein_sequences.fasta" \
    --results_file "$dir_name/temp_work/creating_protein_sequences/fragmented_protein_sequences.fasta"\
    --window_size $window_size --step_size $step_size \
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err" 
elif [ $ROUTE == 3 ]
then 
    # 1. Generate perfonalized protein sequences
        # 1.b running predictions 
    vcf2prot -f $vcf_file -r $fasta_file -vs mt -o "$dir_name/temp_work/creating_protein_sequences/personalized_sequences" \
        >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err" 
    
    # 2. Fragmenting the generated proteins 
    fragment_peptides --input_fasta_dir "$dir_name/temp_work/creating_protein_sequences/personalized_sequences" \
    --results_file "$dir_name/temp_work/creating_protein_sequences/fragmented_protein_sequences.fasta"\
    --window_size $window_size --step_size $step_size --tissue \
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err"
elif [ $ROUTE == 4]
then 
    # 1. Fragmenting the generated proteins 
    fragment_peptides --input_fasta_file $fasta_file \
    --results_file "$dir_name/temp_work/creating_protein_sequences/fragmented_protein_sequences.fasta"\
    --window_size $window_size --step_size $step_size \
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err" 
fi 
## Execute the common Route of each task
#----------------------------------------
# 3. Read the list of alleles and standrdize the output 
allele2standard.py --input_alleles $list_allele --output_results "$dir_name/temp_work/creating_protein_sequences/allele_names_stabdardized.tsv"\
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err" 

# 4. Combine the alleles with the list of fragments
Mixer.py --input_alleles "$dir_name/temp_work/creating_protein_sequences/allele_names_stabdardized.tsv"\
    --input_peptide "$dir_name/temp_work/creating_protein_sequences/fragmented_protein_sequences.fasta" \
    --output_allele "$dir_name/temp_work/creating_protein_sequences/input_to_the_prediction_engine.tsv" \
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err" 

# 5. run the predictions 
prediction_engine_portal.py --base_dir "$dir_name/temp_work/creating_protein_sequences/" \
    --input "$dir_name/temp_work/creating_protein_sequences/fragmented_protein_sequences.fasta" --model_index $model_index \
    --output_path "output/prediction_results.tsv" \
    --unmapped_results "output/unmapped_results.tsv" \
    >> "$dir_name/$( basename $dir_name).logs" 2>> "$dir_name/$( basename $dir_name).err"

# 6. clean up the temp directory
rm -r "$dir_name/temp_work/creating_protein_sequences/" 
rm -r "$dir_name/temp_work/"

echo "INFO:: $(date -u): Exectuion finish" >> "$dir_name/$( basename $dir_name).logs"
# call PIA-inference wrapper engine 
#----------------------------------

