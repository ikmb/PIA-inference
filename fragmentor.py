#!/Users/heshamelabd/miniconda3/envs/ml_env/bin/python3
"""
@brief: A Python tool that takes as input the path to a directory containing FASTA files, it reads all the files, fragment them in parallele using a sliding window approach 
with a step size of 1 and return the generated fragments.
@author: Hesham ElAbd
@contact h.elabd@ikmb.uni-kiel.de
@copyright: Institute of clinical Molecular biology, Univerity of Kiel, Kiel, Germany, 2021  
@Version: 0.0.1
@parallele execution: yes
@parallele mode: Asynchronous multi-processing  
"""
## LOAD THE MODULES
#------------------
import os
import argparse
import sys
import time 
import multiprocessing as mp
import concurrent.futures
from Bio import SeqIO
import pandas as pd 
from typing import Tuple, List, Dict
## DEFINE THE INPUT INTERFACE
#----------------------------
print(f"Execution starts at: {time.ctime()}")
def parse_arguments():
    parser=argparse.ArgumentParser(description=
                            """
                                A Python tool that takes as input the path to a directory containing FASTA files, or a single fasta file. 
                                It reads all the files, fragment them in parallele using a sliding window approach 
                                with a step size of 1 and return the generated fragments.
                            """)
    parser.add_argument('--input_fasta_file', 
                    help = "The path to a FASTA file that can be used for generating input peptides.",
                    type=str, default="?")
    parser.add_argument('--input_fasta_dir',
                    help = "The path to a directory containing multiple FASTA files these files are read and fragmented on parallele.",
                    type=str, default="?")
    parser.add_argument('--window_size',
                    help="The window size to perform the fragmentation",
                     type=int, default=15)
    parser.add_argument('--step_size',
                    help="The step size of the sliding window",
                    type=int, default=1)
    parser.add_argument('--num_worker',
                    help="The number of workers in the worker's pool,defaults to the number of available CPU cores",
                    type=int, default= mp.cpu_count())
    parser.add_argument('--results_path',
                    help="The Path to write the results files as table containing peptide sequences along with the source information.",
                    type=str,default='?') 
    parser.add_argument('--tissue', help="The name of the tissue, only used incase of running predictions against a PIA-M model",   
                    default='-1', type=str, action='store')
    # parse user arguments 
    #---------------------
    args=parser.parse_args()

    ## add args 
    #----------
    if args.input_fasta_file=='?' and args.input_fasta_dir=='?':
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: The input fasta file or input fasta directory has not been provided.\n')
        sys.exit(-1)
    if (not os.path.exists(args.input_fasta_file)) and (not os.path.exists(args.input_fasta_dir)):
        input_name= args.input_fasta_file if args.input_fasta_file!='?' else args.input_fasta_dir
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: The path to the input : {input_name} does not exists.\n')
        sys.exit(-1)
    if args.input_fasta_file!='?' and args.input_fasta_dir!='?':
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: the input fasta file and directory both have been provided. Meanwhile, these variable are mutually exclusive. \n')
        sys.exit(-1)
    # check the value of the window size
    #-----------------------------------  
    if args.window_size not in range(9,22):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: The provided window size is not in range, only values between 9 and 21 are supported.')
        sys.exit(-1)
    
    # Check that we have access write at the resulting directory
    #-----------------------------------------------------------
    if not os.path.exists(os.path.dirname(args.results_path)):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: The base to write the write the results:\
{os.path.dirname(args.results_path)} does not exists\n')
        sys.exit(-1)
    if not os.access(os.path.dirname(args.results_path), os.W_OK):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: you (i.e. user: {os.login()} ) do not have writting access at the results path: {os.path.dirname(args.results_path)}\n')
        sys.exit(-1)

    # return the results
    #-------------------
    return {
        'input_fasta_file':args.input_fasta_file,
        'input_fasta_dir':args.input_fasta_dir,
        'window_size':args.window_size,
        'step_size':args.step_size,
        'results_path':args.results_path,
        'num_worker':args.num_worker,
        'tissue':args.tissue
    }   
   
def fragment(user_argument:Dict[str,str])->pd.DataFrame:
    """A wrapper function used for disbatching the proteome fragmentation function  

    Args:
        user_argument (Dict[str,str]): a dict containing user parameters 

    Returns:
        pd.DataFrame: A pandas dataframe contaning the fragmented proteome results 
    """
    if user_argument['input_fasta_file']!='?':
        fragmented_proteome=fragment_proteome(user_argument['input_fasta_file'],user_argument['window_size'],
            user_argument['step_size'])
    else: 
        fragmented_proteome=fragment_multiple_proteomes(user_argument['input_fasta_dir'],user_argument['window_size'],
            user_argument['step_size'],user_argument['num_worker'])
    if user_argument['tissue']!='-1':
            fragmented_proteome['tissue']=[user_argument['tissue']]*fragmented_proteome.shape[0]
    return fragmented_proteome

def fragment_proteome(path2fasta:str,window_size:int,step_size:int)->pd.DataFrame:
    """Fragment a single proteome and return a table containing all fragmented proteins

    Args:
        path2fasta (str): The path to load the FASTA file 
        window_size (int): The window size 
        step_size (int): The step size of the moving sliding window

    Returns:
        pd.DataFrame: The resulting DataFrame containing peptide sequence 
    """
    # Create the protein and protein dict
    #------------------------------------
    protein_name=path2fasta.split('/')[-1].split('.')[0]
    protein_dict={seq.id.split('|')[1]:str(seq.seq) for seq in SeqIO.parse(path2fasta,'fasta')}
    # Fragment the protein & get the index of each peptide  
    #-----------------------------------------------------
    peptide_sequences=[]
    peptide_names=[]
    for prot_id,prot_seq in protein_dict.items():
        for idx in range(0,len(prot_seq)+1-window_size,step_size):
            peptide_sequences.append(prot_seq[idx:idx+window_size])
            peptide_names.append(protein_name+'$'+prot_id+'_'+str(idx))
    # Prepare the database 
    #---------------------
    results=pd.DataFrame({
        'peptides':peptide_sequences,
        'info':peptide_names
    })
    ## return the results
    #--------------------
    return results

def fragment_multiple_proteomes(path2fasta_dir:str,window_size:int,step_size:int,num_worker:int)->pd.DataFrame:
    """Fragment the proteome of each FASTA in the provided directory on parallel.  

    Args:
        path2fasta_dir (str): The path to a directory containing at least one FASTA file 
        window_size (int): The size of the fragmentation window 
        step_size (int):  The step size of the moving sliding window
        num_worker (int): Number of workers used to run/executing the pipeline 

    Returns:
        pd.DataFrame: The resulting DataFrame containing peptide sequence 
    """
    fasta_files=[os.path.join(path2fasta_dir,fasta_file) for fasta_file in os.listdir(path2fasta_dir) if '.fa' in fasta_file]
    results=pd.DataFrame(columns=['peptides','info','tissue'])
    if fasta_files==[]:
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: The provided input directory: {path2fasta_dir} does not seem to have any Fasta files !!!\n')
        sys.exit(-1)
    else:
        num_worker=min(num_worker,len(fasta_files))
        jobs=[]
        with concurrent.futures.ProcessPoolExecutor(num_worker) as engine:
            for fasta_file in fasta_files:
                jobs.append(engine.submit(
                    fragment_proteome,fasta_file, window_size, step_size
                ) )       
        for fut in concurrent.futures.as_completed(jobs):
            if fut.exception() is not None:
                sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: Fragmenting a proteome gave rise to the following error: {str(fut.exception())}.\n')
                sys.exit(-1)
            proteome_fragments=fut.result()
            results=pd.concat([results,proteome_fragments],axis=0)
    return results
## Create main to load peptides and proteins
if __name__=='__main__':
    user_arguments=parse_arguments()
    proteome_fragments=fragment(user_arguments)
    proteome_fragments.to_csv(user_arguments['results_path'],sep='\t',index=False)