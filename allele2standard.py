#! /usr/bin/env python3
"""
@brief: A python helper tool for standardizing allele names  
@author: Hesham ElAbd
@contact h.elabd@ikmb.uni-kiel.de
@copyright: Institute of clinical Molecular biology, Univerity of Kiel, Kiel, Germany, 2021  
@Version: 0.0.1
@parallele execution: yes
@parallele mode: Asynchronous multi-processing  
"""
## Loading the models
#--------------------
import os
import argparse
import sys
import time 
import pandas as pd 
from typing import * 
#-------------------
## Standardize the allele names 
#------------------------------

def parse_user_arguments()->Dict[str,str]:
    """ parse user input 

    Returns:
        Dict[str,str]: the parsed user arguments 
    """
    parser=argparse.ArgumentParser(description="Correcting allele names into a form supported or understood by the prediction engine")
    parser.add_argument('--input_table', help="The input list of alleles, a table with one columns", default='?',type=str)
    parser.add_argument('--output_table', help="A table with one columns representing the correct allele names", default='?',type=str)
    # parse the argument and check the validity of the the input 
    args=parser.parse_args()
    # check that the input path is valid 
    if args.input_table=='?':
        sys.stderr.write(f'{time.ctime()} CRITICAL:: The base directory has not been provided. Exiting\n')
        sys.exit(-1)
    # check that the output path is valid 
    if not os.path.exists(os.path.dirname(args.output_table)):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: The base to write the write the results:\
{os.path.dirname(args.output_table)} does not exists\n')
        sys.exit(-1)
    if not os.access(os.path.dirname(args.output_table), os.W_OK):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: you (i.e. user: {os.login()} ) do not have writting access at the results path: {os.path.dirname(args.results_path)}\n')
        sys.exit(-1)
    # Return the results
    #-------------------
    return {
        'input_table':args.input_table,
        'output_table':args.output_table
    }

def correct_allele_name(input_path:str)->pd.DataFrame:
    """Correct the allele names

    Args:
        input_path (str): The path to the input list of alleles.

    Returns:
        pd.DataFrame: a table containing corrected list of alleles.
    """
    input_table=pd.read_csv(input_path,sep='\t')
    new_allele_list=[]
    for allele in input_table.iloc[:,0].to_list():
        if '/' not in allele:
            sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: can not standardize input allele names {allele}\n')
            sys.exit(-1)
        if 'DR' in allele:
            new_allele_list.append(allele.split('/')[1].strip('HLA-').replace('*','_').replace(':',''))
        elif 'DP' in allele or 'DQ' in allele:
            alpha,beta=allele.split('/')
            alpha=alpha.replace('*','').replace(':','')
            beta=beta.strip('HLA-').replace('*','').replace(':','')
            new_allele_list.append('alpha'+'-'+'beta')
    # Return a new table to hold the results
    new_table=pd.DataFrame({'alleles':new_allele_list})
    # return the table 
    return new_table

## Executing the main logic of the pipeline 
#------------------------------------------
if __name__=='__main__':
    user_arguments=parse_user_arguments()
    new_table=correct_allele_name(user_arguments['input_table'])
    new_table.to_csv(user_arguments['output_table'],index=False, sep='\t')
    

    