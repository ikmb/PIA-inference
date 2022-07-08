#!/usr/bin/env python3
"""
@author: Hesham ElAbd
@brief: Create protein sequences a table containing protein sequence
"""
## Load modules
#--------------
import os 
import sys
import argparse
import time
import pandas as pd 
from Bio import SeqIO
from tqdm import tqdm
from typing import * 
## define the databases
#----------------------
def parse_argument()->Dict[str,Any]:
    """A helper function to parse user inputs 

    Returns:
        Dict[str,Any]: the argument name and value
    """
    parser=argparse.ArgumentParser(description='Mut2Prot, generates a table of protein sequences using a reference input sequence and a collection of mutation')
    
    parser.add_argument('--input_table',help='A table of three columns, the first columns represents the uniprot id, the second represents the genetic mutation, \
        expressed as amino_acid_position_amino_acid, e.g. K341L or N211T, etc. The position **MUST BE ONE INDEXED**. The 3rd columns contains the phenotypes (Optional). Incase you genetic data is more complex than\
            this representation. Please reformat your data as a VCF files and use VCF2Prot for execution. **FILES MUST BE TAB SEPERATED FILES**',
            default='?',action='store')

    parser.add_argument('--input_fasta', help="The path to a fasta file where the header file is ONLY composite of the\
         FASTA sequence ID. Also, the ID used in the gentic table matches the ID used in the gentic table.\
         Otherwise a mismatch must be conducted, followed by the reference protein sequence.",
         default='?',action='store')
    
    parser.add_argument('--window_size', help="The window size for fragmentation, can be an integer bigger than or equal to 9 and smaller than or equal to 21",
        default=15,type=int,action='store')
    
    parser.add_argument('--results_path', help="The path to write the generated protein sequence to a FASTA file.",
        default='?', action='store'
    )
    # Parse the user arguments
    args=parser.parse_args()
    # Check the validity of the input table
    #-------------------------------------- 
    if args.input_table=='?':
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The input genetic table has not been provided.\n')
        sys.exit(-1)
    if not os.path.exists(args.input_table):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The path to the input table: {args.input_table} has not been provided.')
        sys.exit(-1)
    try:
        dev_read=pd.read_csv(args.input_table,nrows=10,sep='\t')
    except Exception as exp:
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: Reading the input table provided at: {args.input_table}\
             risen the following exception: {str(exp)}')
        sys.exit(-1)
    if dev_read.shape[1] not in [2,3]:
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: Incorrect table provided at: {args.input_table}\
             expect number of columns to 2 or 3 columns. However, your input has : {dev_read.shape[1]} columns')
        sys.exit(-1)
    # Check the validity of the resulting input FASTA files
    #-------------------------------------------------------
    if args.input_fasta == '?':
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The input FASTA has not been provided.\n')
        sys.exit(-1)
    if not os.path.exists(args.input_table):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The path to the input FASTA file: {args.input_table} has not been provided.')
        sys.exit(-1)
    # Check that we have access write at the resulting directory
    #-----------------------------------------------------------
    if not os.path.exists(os.path.dirname(args.results_path)):
        print(f"Base dir is: {os.path.dirname(args.results_path)}, Does it exist: {os.path.exists(os.path.dirname(args.results_path))}")
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The base to write the write the results:\
             {os.path.dirname(args.results_path)} does not exists.')
        sys.exit(-1)
    if not os.access(os.path.dirname(args.results_path), os.W_OK):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: you (i.e. user: {os.login()} ) do not have writing access at the results path: {os.path.dirname(args.results_path)}\n')
        sys.exit(-1)
    # check the value of the window size
    #-----------------------------------  
    if args.window_size not in range(9,22):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The provided window size is not in range, only values between 9 and 21 are supported.')
        sys.exit(-1)
    ## return the results
    #--------------------
    return {
        'genetic_table':args.input_table,
        'input_fasta':args.input_fasta,
        'results_path':args.results_path,
        'window_size':args.window_size,
    }

def load_fasta_sequence(path2fasta:str)->Dict[str,str]:
    """Load FASTA sequences

    Args:
        path2fasta (str): The path to load the input file

    Returns:
        Dict[str,str]: a dict of sequence id to sequence name 
    """
    results=dict()
    for seq in SeqIO.parse(path2fasta,'fasta'):
        results[seq.id]=str(seq.seq)
    return results

def generate_personalized_proteins(proteome:Dict[str,str],path2gentic_table:str,window_size:int)->pd.DataFrame:
    """Load personalized proteins

    Args:
        proteome (Dict[str,str]): A table contain protein names as keys and protein sequence as values. 
        path2gentic_table (str): The path to genetic table.

    Returns:
        pd.DataFrame: The generated data frame. 
    """
    # Load genetic table
    #-------------------
    genetic_table=pd.read_csv(path2gentic_table,sep='\t')
    if genetic_table.shape[1]==2:
        genetic_table['Pheno']=['UNK']*genetic_table.shape[0]
    # allocate list to hold the results
    #----------------------------------
    generated_results=pd.DataFrame(columns=['peptide_sequence','peptide_index','Pheno'])
    # loop over all the elements in the input table
    #----------------------------------------------
    for row in tqdm(genetic_table.itertuples()):
        try: 
            protein_seq=proteome[row.name]
        except Exception as exp:
            sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encountered the following problem: The sequence of protein {row.name} is not defined the input proteome.')
            sys.exit(1)
        # obtain the peptides and name
        mutated_peptides=apply_genetic_variant(row.genetic_variant, protein_seq, row.name, window_size)
        mutated_peptides['Pheno']=[row.Pheno]*mutated_peptides.shape[0]
        # append the results along the x-axis
        generated_results=pd.concat([generated_results,mutated_peptides],axis=0)
    # Fragment the reference proteome 
    #--------------------------------
    reference_proteome_peptides=fragment_reference_proteome(proteome,window_size)
    # combine the results of the two datasets, i.e. the reference and the mutations
    #------------------------------------------------------------------------------
    generated_results=pd.concat([generated_results,reference_proteome_peptides],axis=0)
    # return the results 
    return generated_results

def fragment_reference_proteome(input_proteome:Dict[str,str], window_size:int)->pd.DataFrame:
    """Fragment the reference proteomes into overlapping peptides 

    Args:
        input_proteome (Dict[str,str]): The input proteome where keys are protein names while values are protein sequences. 
        window_size (int): The size of the enerated fragment. 

    Returns:
        pd.DataFrame: The resulting fragmented proteins. 
    """
    peptide_seqs=[]
    peptide_index=[]
    Pheno=[]
    for protein_name, protein_sequences in input_proteome.items():
        for idx in range(0,len(protein_sequences)+1-window_size,1):
            peptide_seqs.append(protein_sequences[idx:idx+window_size])
            peptide_index.append(protein_name+f'_Reference_'+str(idx))
            Pheno.append('Reference')
    # Make a dataframe of the generated fragments and return it 
    return pd.DataFrame({'peptide_sequence':peptide_seqs,'peptide_index':peptide_index,'Pheno':Pheno})

def apply_genetic_variant(genetic_mutation:str, protein_sequence:str, protein_name:str, window_size:int)->pd.DataFrame:
    """Apply genetic mutations to the input protein sequence and return the results

    Args:
        genetic_mutation (str): The genetic name of the mutation, e.g. K102L 
        protein_sequence (str): The sequence of the protein
        protein_name (str): The name of the protein 
        window_size (int): The window size, e.g 15 amino acid

    Returns:
        pd.DataFrame: Return a data frame containing all generated peptide fragments along with a name and an index
    """
    try:
        _, position, altered_aa = parse_genetic_code(genetic_mutation)
    except Exception as exp:
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: parsing the following genetic code: {genetic_mutation} gave rise to the following exception: {str(exp)}.')
        sys.exit(1)
    if position > len(protein_sequence):
        sys.stderr.write(f'ERROR: {time.ctime()}: I {sys.argv[0]} encounterred the following problem: the position of the mutation: {position} is longer than protein: {protein_name} length: {len(protein_sequence)}.')
        sys.exit(1)
    # extract the upper and lower bound for fragmentation
    #----------------------------------------------------
    start, end= max(0,position+1-window_size),min(len(protein_sequence),position+window_size)
    # Generate the altered protein sequence
    #--------------------------------------
    altered_protein=protein_sequence[:position]+altered_aa+protein_sequence[position+1:]
    target_sequence=altered_protein[start:end]
    # loop over all alleles 
    #----------------------
    peptide_seqs=[]
    peptide_index=[]
    for idx in range(0,len(target_sequence)+1-window_size,1):
        peptide_seqs.append(target_sequence[idx:idx+window_size])
        peptide_index.append(protein_name+f'_{genetic_mutation}_'+str(idx))
    # Create a data frame of the results
    #-----------------------------------
    results_df=pd.DataFrame({'peptide_sequence':peptide_seqs,'peptide_index':peptide_index})
    ## Return the results data frame 
    #-------------------------------
    return results_df

def parse_genetic_code(genetic_mutation:str)->Tuple[str,int,str]:
    """parse the genetic mutation to extract its relevent information content. 

    Args:
        genetic_mutation (str): The genetic mutation string

    Returns:
        Tuple[str,int,str]: a tuple of three elements containing the reference amino acids, the position and the altered amino acids
    """
    reference=genetic_mutation[0]
    altered=genetic_mutation[-1]
    position=int(genetic_mutation[1:-1])-1
    return reference, position, altered

## Start writing the main execution logic of the script
if __name__=='__main__':
    user_argument=parse_argument()
    protein_sequences=load_fasta_sequence(user_argument['input_fasta'])
    generated_personalized_sequences=generate_personalized_proteins(protein_sequences,user_argument['genetic_table'],user_argument['window_size'])
    generated_personalized_sequences.to_csv(user_argument['results_path'],sep='\t',index=False)
