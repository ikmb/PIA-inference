#! /usr/bin/env python3
"""
@author: Hesham ElAbd
@brief: Mixing the list of alleles with the list of input peptides 
"""
## Load the modules
#------------------
import time 
import os
import sys
import argparse
import pandas as pd
from tqdm import tqdm 
from typing import * 
## Load the modules
#------------------
def parse_user_inputs()->Dict[str,str]:
    """Parse use input and returns a dict containing user parameters. 

    Returns:
        Dict[str,str]: a dict containing user provided elements
    """
    # Define the parse to parse the results
    #--------------------------------------
    parser=argparse.ArgumentParser(description='Mixer, generates a table of containing a the product of multiplying the set of peptides and the set of supported alleles')
    
    parser.add_argument('--input_peptides',help='A table of two columns containing the peptide sequence and the allele name',
            default='?', action='store')

    parser.add_argument('--standardized_allele_name',help="Standardized allele name, a table with one columns containing a header along with a collection of allele names",
    default='?', action='store')

    parser.add_argument('--results_path',help="The output path to write the results table.",
     default='?', action='store'
    )
    # Parse user arguments
    #---------------------
    args=parser.parse_args()
    #--------------
    # Load input genetic table and check for correctness 
    if args.input_peptides=='?':
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: The input list of peptides has not been provided.\n")
        sys.exit(-1)
    try:
        dev_read=pd.read_csv(args.input_peptides,nrows=10,sep='\t')
    except Exception as exp:
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: Reading the input genetic table provided at: {args.input_peptides}\
             risen the following exception: {str(exp)}")
        sys.exit(-1)
    if dev_read.shape[1] != 2:
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: Incorrect table provided at: {args.input_peptides}\
             expect number of columns is 2. However, your input has : {dev_read.shape[1]} columns")
        sys.exit(-1)
    
    # Load the list of alleles
    #-------------------------
    if args.standardized_allele_name=='?':
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: The list of alleles has not been provided.\n")
        sys.exit(-1)
    if not os.path.exists(args.standardized_allele_name):
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: The path to the list of alleles: {args.input_peptides} has not been provided.")
        sys.exit(-1)
    ## check the output is path is valid
    #-----------------------------------
    if not os.path.exists(os.path.dirname(args.results_path)):
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: The base to write the results:\
             {os.path.dirname(args.results_path)} does not exists.")
        sys.exit(-1)
    if not os.access(os.path.dirname(args.results_path), os.W_OK):
        sys.stderr.write(f"ERROR: {time.ctime()}: I {sys.argv[0].split('/')[-1]} encountered the following problem: you (i.e. user: {os.login()} ) do not have writting access at the results path: {os.path.dirname(args.results_path)}\n")
        sys.exit(-1)
    ## return the results
    #--------------------
    return {
        'input_peptides':args.input_peptides,
        'standardized_allele_name':args.standardized_allele_name,
        'results_path':args.results_path
    }

def mix_inputs(path2pep:str, path2alleles:str)->pd.DataFrame:
    """ Loads the input table of peptides along with the list of alleles and mix them together to generate the input table for prediction

    Args:
        path2pep (str): The path to load the input peptide table 
        path2alleles (str): The path to load the path to alleles

    Returns:
        pd.DataFrame: The resulting dataframe containing a combination of all peptides against all alleles
    """
    ## Load the two table 
    #--------------------
    peptides=pd.read_csv(path2pep,sep='\t')
    alleles=pd.read_csv(path2alleles,sep='\t')
    # allocate the results table
    result_df=pd.DataFrame(columns=list(peptides.columns)+['allele'])
    # repeat the table for each input allele
    for allele in alleles.iloc[:,0].to_list():
        temp_ds=peptides.copy(deep=True)
        temp_ds['allele']=[allele]*temp_ds.shape[0]
        # append the results to the table 
        result_df=pd.concat([result_df,temp_ds],axis=0)
    # return the results 
    #-------------------
    return result_df

## add the main execution logic 
if __name__=='__main__':
    user_arguments=parse_user_inputs()
    input2predictions=mix_inputs(user_arguments['input_peptides'],user_arguments['standardized_allele_name'])
    input2predictions.to_csv(user_arguments['results_path'],index=False,sep='\t')
