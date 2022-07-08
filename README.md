# Peptide immune annotator pipeline

## Aims

To provides a framework for modeling peptide HLA-II at large-scale using different input types.

## Implementation

The pipeline is currently implemented as a Bash script that first, glues together different Python scripts, second parse user argument and select the correct execution path. Finally, it manages the logging and the correct execution of each script.

### Building blocks

The pipeline is built using different modular units, these are:

#### I. *Framgentor*

The *fragmentor* is a Python-based command line tool (CLT) that is used for fragmenting input proteins into peptides with a fixed size using a sliding window approach. These can either be a single FASTA File or multiple FASTA files. Incase of a single file, execution is single-threaded. Meanwhile, with a list of FASTA files a pool of worker threads is used for parsing the input. Regardless, of the execution methods, peptides generate from a single or multiple proteomes are written as a single tsv file. An example for using the fragmentor as a standalone tool is show below

```bash
# 1. make sure you have either cloned the repository with 
git clone https://github.com/ikmb/PIA-inference

# 2. Change the directory into the PIA-interface directory 
cd PIA-inference

# 3. create the conda environment and install dependencies as described below, if done already please jump to the next step 

# 4. let's call the fragmentor script to see how it is working 
./fragmentor.py -h 

# 5. Then let's call with the proteome of the SARS-COV 2 which was download from Uniprot database 
./fragmentor.py --input_fasta_file examples/SARS_COV2.fasta \
    --window_size 15 \
    --step_size 1 \
    --num_worker 1 \
    --results_path examples/fragmentation_results.tsv 

# 6. This should take a few seconds to finish (depending on your machine) but after it finished we can check for the file as follow
head examples/fragmentation_results.tsv 

# 7. If we have a collection of files we create a directory to hold them 

mkdir examples/example_proteins

# 8. We copy all fasta file to the example proteins 
cp examples/*fasta example_proteins

# 9. Call the program  
./fragmentor.py --input_fasta_dir examples/SARS_COV2.fasta \ # Notice the change in parameter from input_fasta_file to input_fasta_dir
    --window_size 15 \
    --step_size 1 \
    --num_worker 1 \
    --results_path examples/fragmentation_results_all.tsv 

# 10. look at the first few lines of the generated files 
head examples/fragmentation_results_all.tsv 

```

#### II. *Mut2Prot*

It is a Python-bashed CLT that is parsing a tsv table containing genetic mutations in order to generate a personalized or sample-specific proteomes, *i.e.* a collection of protein sequences.

```bash
# 1. make sure you have either cloned the repository with 
git clone https://github.com/ikmb/PIA-inference

# 2. Change the directory into the PIA-interface directory 
cd PIA-inference

# 3. create the conda environment and install dependencies as described below, if done already please jump to the next step 

# 4. let's call Mut2prot to generate the list of peptides 

./Mut2Prot.py --input_table examples/input_genetic_variant_corrected.tsv\ # the input table 
            --input_fasta examples/example_protein.fasta\ # the input fasta file containing reference proteins
            --window_size 15\ # the fragmentation window 
            --results_path examples/results_after_collecting_genetic_data.tsv # the generated results 
```

#### III. *vcf2prot*

VCF2Prot is a Rust based executable that can be used for the large scale generation of sample-specific proteomes from VCF files, for more information regarding the use of VCF2prot please [click here](https://github.com/ikmb/vcf2prot/tree/main)

#### IV. *allele2standard*

It is a Python-based command line tool that takes a list of alleles written in the standard notation as an input and return formatted names that is used for handling extracting pseudo-sequences. An example of using the pipeline is sen below

```bash

# 1. let's write a group of alleles to a file to prepare an input case
echo "List_alleles" >> examples/test_case_allele2standard.tsv  
echo  "HLA-DRA1*01:01/HLA-DRB1*15:01"  >> examples/test_case_allele2standard.tsv 
echo  "HLA-DRA1*01:01/HLA-DRB1*01:01"  >> examples/test_case_allele2standard.tsv  
echo  "HLA-DRA1*01:01/HLA-DRB1*13:01"  >> examples/test_case_allele2standard.tsv  

# 2. let's call the program 
./allele2standard.py --input_table examples/test_case_allele2standard.tsv\
                     --output_table examples/test_case_allele2standard_output.tsv 

# 3. let's read the output of the file 
cat examples/test_case_allele2standard_output.tsv

# this shall print the following 4 four lines 
alleles
DRB1_1501
DRB1_0101
DRB1_1301
```

### Mixer

The Mixer is a CML line that is used for mixing the list of alleles with the input peptides to prepare the input to the prediction model as explained below

### The prediction engine

This where the interaction between the generated peptides and HLA-II proteins is estimated using PIA-S and PIA-M, a short description for using the prediction engine is illustrated below.
Both relies heavy on the OmLiT library. [Clink here](https://github.com/ikmb/OmLiT) to read more about OmLiT and to install it on your local machine.

## List of dependencies

### I. Python dependencies

#### A. Install with Pip

```bash
# we highly recommend to create a new environment using Conda to have a clean and isolated environment 
conda create -n pia_inf -y # here we chose pia_inf but you can replace it with any name you like 

# Activate the Conda environment 
conda activate pia_inf

# Install pip 
conda install -y pip 

# install the dependencies 

pip install -r python_requirements.txt 

```

#### B. Install with *Conda*

```bash
# we highly recommend to create a new environment using Conda to have a clean and isolated environment 
conda create -n pia_inf -y --fille python_requirements_conda.txt # here we chose pia_inf but you can replace it with any name you like 

# Activate the Conda environment 
conda activate pia_inf
```

## Using the list of the pipeline

### I. calling the pipeline for help messages

```bash
./PIA_inference.sh -h
```

### II. Testing the first execution mode

Here, we assume that we have a table with input peptides and list of alleles to be fed directly to the prediction engine as described below

```bash
# 1. call the pipeline 
./PIA_inference.sh -d tests/test_one/ -s test_input_table.tsv -m 1 -w 15 -z 1

# 2. Incase any error was encountered the output will be directed to the errors file, we cat read it as follow 
cat tests/test_one/test_one.err 
# or using 
cat tests/test_one/test_one.err  | less 

# 3. we can also look at any output message printed by the pipeline using
cat tests/test_one/test_one.logs

# 4. finally, the output file will be located at: 

head tests/test_one/output/prediction_results.tsv
```

### III. Testing with the second execution model

Here, the aim is to run predictions using a file containing FASTA files and a file containing simple genetic mutations against a list of alleles

```bash

## NOTE 
#------
# Skip the following point if you have already run it before with testing allele2standard

# 1. let's prepare the list of alleles 
echo "List_alleles" >> examples/test_case_allele2standard.tsv  
echo  "HLA-DRA1*01:01/HLA-DRB1*15:01"  >> examples/test_case_allele2standard.tsv 
echo  "HLA-DRA1*01:01/HLA-DRB1*01:01"  >> examples/test_case_allele2standard.tsv  
echo  "HLA-DRA1*01:01/HLA-DRB1*13:01"  >> examples/test_case_allele2standard.tsv  

# 2. let's call the pipeline 
./PIA_inference.sh -d tests/test_two -g input_genetic_variant_corrected.tsv \
                -f example_protein.fasta \
                -a test_case_allele2standard.tsv \
                -m 1 -w 15 -z 1
```

### IV. Testing with the 3rd execution model

Here, the aim is to run the predictions using a VCF file containing patient genetic variants, a reference proteome stored in FASTA file and a list of alleles stored in a text file.

#### Notes

1. Please make sure vcf2prot is installed and working on your system before continuing, this can be done using:  

```bash
./vcf2prot -h 
```

for more information regarding the installation check the webpage above, here we are going to utilize a VCF file containing genetic mutations observed across sample HG00096. This dataset was obtained from the [1000Genome project](http://www.internationalgenome.org). While the reference proteome will be obtained from Ensemble reference's specifically, the protein sequence of each transcript.

#### Calling the pipeline

```bash

# 1. let's first decompress the generated VCF file 
gunzip vcf_dev.vcf.gz 

# 2. let's de compress the reference sequence 
gunzip reference_proteome.fasta.gz

# 3. call the program 
./PIA_inference.sh -d tests/test_three/ \
                -f reference_proteome.fasta \
                -a test_case_allele2standard.tsv \
                -v vcf_dev.vcf -m 1 -w 15 -z 1
```

### IV. Testing with the 4th execution model

Here, the aim is to execute the model with a FASTA file representing a proteome of interest, here we are going to focus on COVID-19 proteome as we did before

```bash
# 1. we can just call the program 

./PIA_inference.sh -d tests/test_four/ \
    -f SARS_COV2.fasta \
    -a test_case_allele2standard.tsv\
    -m 1 -w 15 -z 1
```

## Funding

The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’)
