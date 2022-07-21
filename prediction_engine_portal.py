#!/usr/bin/env python3
"""
@author: Hesham ElAbd
@brief: Predicting the interaction between an input list of peptides and a list of HLA proteins.    
@version: 0.1.0 pre-alpha. 
@execution: Heterogenous using CPUs and GPUs.
@preprocessing: OmLiT for sequencing encoding an omics preprocessing. 
@copyright: institute of clinical molecular biology (IKMB), Kiel, Germany. 
@data: 12.03.2022
"""
## import the modules 
#--------------------
import pickle
import os
import sys
import argparse
import time 
import pandas as pd
import numpy as np 
import OmLiT as linker 
import tensorflow as tf
from ImFormers import * 
from typing import * 
import logging
from tqdm import tqdm 
from Bio import SeqIO
## Set the logger to a show only errors
#--------------------------------------
tf.get_logger().setLevel(logging.ERROR)
## define the analysis function 
#--------------------------
def parse_user_inputs()->Dict[str,str]:
    """Parse use input and returns a dict containing user parameters 

    Returns:
        Dict[str,str]: a dict containing user provided user parameters  
    """
    # Define the parse to parse the results
    #--------------------------------------
    parser=argparse.ArgumentParser(description="Predicting the interaction between an input list of peptides and a list of HLA proteins using either PIA-S or PIA-M")

    parser.add_argument('--base_dir',help="The base directory for running the prediction job, it contains the user input and output directories.",
        default='?',type=str,action='store'
    )
    parser.add_argument('--input',help="The input table containing three columns, the first is the peptide sequence, the second is the allele and the third is the tissue name.\
    The first two columns are mandatory inputs while the third is optional, incase PIA-S is used then the the first two column are not used meanwhile, if PIA-M is used then the default expression of total PBMCs is used.\
        The input table must be must be TAB seperated.",
    default='?',type=str,action='store'
    )
    parser.add_argument('--model_index', 
    help="The index of the model, 1. is PIA-S trained on public databases, 2. PIA-S trained on public and inhouse datasets, 3. is PIA-M trained on public databases and 4.\
    is PIA-M trained on all public and inhouse databases.",
    default=1,type=int,action='store'
    )
    parser.add_argument('--output_path',
    help="The path to write the output results, defaults to prediction_results.tsv in the current working directory",
    default=os.path.join(os.getcwd(),'output/prediction_results.tsv'),
    type=str
    )
    parser.add_argument('--unmapped_results',
    help="The path to write the unmapped results, defaults to unmapped_results.tsv in the current working directory.",
    default=os.path.join(os.getcwd(),'output/unmapped_results.tsv'),
    type=str
    )
    parser.add_argument('--tissue', 
    help="The name where peptide presentation to be studied",
    default=-1,
    type=int
    )
    ## Parse input parameters
    #-----------------------
    args=parser.parse_args()
    ## check that provided user parameters and update the output
    #-----------------------------------------------------------
    if args.base_dir=='?':
        sys.stderr.write(f'{time.ctime()} CRITICAL:: The base directory has not been provided. Exiting\n')
        sys.exit(-1)
    elif not os.path.exists(args.base_dir):
        sys.stderr.write(f'{time.ctime()} CRITICAL:: The provided base directory does not exist. Exiting \n')
        sys.exit(-1)
    else:
        base_dir=args.base_dir
    ## check that we have sub-folder to hold the output 
    try: 
        os.mkdir(os.path.join(base_dir,'output'))
    except FileExistsError:
        pass
    # Load input table 
    if args.input=='?': 
        sys.stderr.write(f'{time.ctime()} CRITICAL:: The input file has not been provided.\n')
        sys.exit(1)
    elif not os.path.exists(args.input):
        sys.stderr.write(f'{time.ctime()} CRITICAL:: The provided input file does not exist!!\n')
        sys.exit(1)
    else: 
        input_path=args.input
    # check the model index
    if args.model_index not in  [1,2]:
        sys.stderr.write(f"{time.ctime()} CRITICAL:: The model index is not valid, the following four models are only valid: 1. is PIA-S, 2. PIA-M. While you input is: {args.model_index} \n")
        sys.exit(1)
    else: 
        model_index=args.model_index
    # Check the validity of the output result file 
    if os.access(os.path.dirname(args.output_path), os.W_OK):
        if base_dir=='.':
            results_path=args.output_path
        else:
           # print(f"The base dir is: {base_dir}")
            results_path=os.path.join(os.path.dirname(args.output_path),base_dir.strip('/').split('/')[-1]+'.prediction_results.tsv')
    else:
        print(f"The provide path is: {args.output_path} and the base directory is: {os.path.dirname(args.output_path)}")
        sys.stderr.write(f"{time.ctime()} CRITICAL:: invalid writhing path to write the output results, user: {os.getlogin()} does not not have writhing access at: {args.output_path}, are you sure the output path is valid. \n")
        sys.exit(-1) 
    # Check the validity of the output result file 
    if os.access(os.path.dirname(args.unmapped_results), os.W_OK):
        if base_dir==".":
            unmapped_results=args.unmapped_results
        else:
            unmapped_results=os.path.join(os.path.dirname(args.unmapped_results),base_dir.strip('/').split('/')[-1]+'.unmapped_results.tsv')
    else:
        sys.stderr.write(f"{time.ctime()} CRITICAL:: invalid writing path to write the unmapped results, user: {os.getlogin()} does not not have writting access at: {args.unmapped_results}, are you sure the output path is valid. \n")
        sys.exit(-1) 

    # return the results
    #-------------------
    return {
        'base_dir':base_dir,
        'input':input_path,
        'model_index':model_index,
        'output_path':results_path,
        'unmapped_results':unmapped_results,
        'tissue':args.tissue
    }

def load_list_of_tissues()->List[str]:
    """Load the list of tissues from the assets

    Returns:
        List[str]: the list of supported tissue names 
    """
    return pd.read_csv('/work_ifs/sukmb418/PIA-inference/assets/list_unique_tissue.txt',sep='\t',header=None).iloc[:,0].to_list()

def get_tissue_name(tissue_index:int,list_of_tissues:List[str])->str:
    """Get the tissue name from the list of tissues 

    Args:
        tissue_index (int): The index of the tissue 
        list_of_tissues (List[str]): the list of supported tissues

    Returns:
        str: the tissue name
    """
    # Get the index of the tissue 
    #----------------------------
    if tissue_index-1>len(list_of_tissues): 
        sys.stderr.write(f"{time.ctime()} ERROR: invalid index, you index is {tissue_index} for a list of size {len(list_of_tissues)}")
        sys.exit(-1)
    elif tissue_index==0:
        sys.stderr.write(f"{time.ctime()} ERROR: invalid index, you index is {tissue_index} for a list that is zero-indexed")
        sys.exit(-1)
    elif tissue_index==-1:
        return 'total PBMC'
    else: 
        return list_of_tissues[tissue_index-1]

def load_model(model_index:int)->tf.keras.models.Model:
     """Load and parse the model and other assets

     Args:
         model_index (int): The index of the model: 1. is PIA-S trained on public databases, 2. PIA-S trained on public and inhouse datasets, 3. is PIA-M trained on public databases and 4.\
     is PIA-M trained on all public and inhouse databases. 

     Returns:
         tf.keras.models.Model: A tensorflow.keras.model that match the load index
     """
     mirrored_strategy = tf.distribute.MirroredStrategy()
     with mirrored_strategy.scope():
          if model_index==1:
               model=create_basic_imformer(maxlen=55,vocab_size=28,embedding_dim=32,num_heads=4,feedforward=64,num_blocks=3)
               model.load_weights('/work_ifs/sukmb418/PIA-inference/assets/PIA_S.h5')
               print(f"{time.ctime()} INFO:: PIA-S has been loaded")
          elif model_index==2:
               model=create_subcellular_location_transcription_contexted_dist_to_gly_imformer(maxlen=55,vocab_size=28,embedding_dim=32,num_heads=4,feedforward=64,num_blocks=3) 
               model.load_weights('/work_ifs/sukmb418/PIA-inference/assets/PIA_M.h5')
               print(f"{time.ctime()} INFO:: PIA-M has been loaded")
          else: 
               raise ValueError(f"{time.ctime()}:: Unsupported model index, current model only support 4 indices,\
                    1 --> PIA-S (public), 2 --> PIA-S (public+inhouse), 3 --> PIA-M (public), 4 --> PIA-M (public+inhouse). while your input is: {model_index}")
     return model

def encode_input_for_pia_s(path2input:str)->Tuple[pd.DataFrame,np.ndarray,pd.DataFrame]:
    """prepare and encode the data for running inferences using PIA-S

    Args:
        path2input (str): The path to load the input file 
    Raises:
        IOError: incase reading the file failed (1) or if it contain in correct number of alleles

    Returns:
        Tuple[pd.DataFrame,np.ndarray,pd.DataFrame]: A tuple of three containing
        1. Pandas dataframe containing the mapped peptides, composite of two columns, the peptide and the HLA allele.
        2. np.ndarray with the following shape: (num_examples,55) --> input to PIA
        3. pandas dataframe containing the unmapped inputs, containing peptide, allele and the tissue. 
    """
    ### Load the input table of peptides and proteins     
    #------------------------------------------------
    try:
        input_table=pd.read_csv(path2input,sep='\t')
        if 'info' in input_table.columns and input_table.shape[1]==3:
            input_table.columns=['peptide','info','allele']
        elif 'peptide_index' in input_table.columns and input_table.shape[1]==3:
            input_table.columns=['peptide','info','allele']
        elif input_table.shape[1]==3:
            input_table.columns=['peptide','allele','tissue']
        elif input_table.shape[1]==2:
            input_table.columns=['peptide','allele']
        elif input_table.shape[1]==4:
            input_table.columns=['peptide',	'info','phenotype', 'allele']
        else:
            raise ValueError(f"In correct number of input columns, excepted 2 or 3 columns, however, the input has: {input_table.shape[0]}")
    except Exception as exp:
        raise IOError(f"{time.ctime()} CRITICAL:: Reading the input table failed with the following error: {exp}")
    
    ### Load the pseudo sequence
    #---------------------------
    pseudo_sequences=pd.read_csv('/work_ifs/sukmb418/PIA-inference/assets/pseudosequence.2016.all.X.dat',header=None,sep='\t')
    pseudo_sequences.columns=['Allele','pseudo_sequence']

    ### get unmapped due-to un matached alleles
    #------------------------------------------
    pseudo_sequences_look_up=pd.read_csv('/work_ifs/sukmb418/PIA-inference/assets/pseudosequence.2016.all.X.dat',header=None,sep='\t')
    pseudo_sequences_look_up.columns=['Allele','pseudo_sequence']
    ### Remove un filtered 
    print(f"{time.ctime()} INFO:: filtering the input table with {input_table.shape[0]} elements.")
    unmapped=input_table.loc[~((input_table.peptide.str.isalpha()) & (input_table.allele.isin(pseudo_sequences_look_up.Allele))),]
    mapped=input_table.loc[((input_table.peptide.str.isalpha()) & (input_table.allele.isin(pseudo_sequences_look_up.Allele))),]
    print(f"{time.ctime()} INFO:: finished with filtering the inputs, number of unmapped is: {unmapped.shape[0]} while number of mapped is: {mapped.shape[0]}")
    
    ## Extract the peptide and the HLA
    #---------------------------------
    peptides =mapped.peptide.to_list()
    HLA=mapped.allele.to_list()
    ### Encoding input peptides:
    #---------------------------
    print(f"{time.ctime()} INFO:: Extracted the pseudo-sequences for each allele")
    #-------------------
    encoded_input,_=linker.annotate_and_encode_input_for_pia_s((peptides,HLA), 21, '/work_ifs/sukmb418/PIA-inference/assets/pseudosequence.2016.all.X.dat')

    ### Make a dataframe from the results
    #------------------------------------
    if 'tissue' in input_table.columns and 'phenotype' not in input_table.columns :
        # Build a dict to lookup the tissue of the mapped peptides
        #---------------------------------------------------------
        print(f"INFO: {time.ctime()} :: Building a look up table between peptides and tissues")
        tissue_lookup={peptide:tissue for peptide,tissue in tqdm(zip(input_table.peptide,input_table.tissue))}
        # Build a table containing the results
        #----------------
        input_peptide=pd.DataFrame({
            'peptide':peptides,
            'allele':HLA,
            'tissue':[tissue_lookup[pep] for pep in peptides]
    })
    elif 'phenotype' in input_table.columns and 'tissue' not in input_table.columns : 
        # Build a dict to lookup the tissue of the mapped peptides
        #---------------------------------------------------------
        print(f"INFO: {time.ctime()} :: Building a look up table between peptides and tissues")
        phenotype_lookup={peptide:phenotype for peptide,phenotype in tqdm(zip(input_table.peptide,input_table.phenotype))}
        print(f"INFO: {time.ctime()} :: Building a look up table between peptides and peptide information")
        info_lookup={peptide:phenotype for peptide,phenotype in tqdm(zip(input_table.peptide,input_table['info']))}
        # Build a table containing the results
        #----------------
        input_peptide=pd.DataFrame({
            'peptide':peptides,
            'allele':HLA,
            'Phenotype':[phenotype_lookup[pep] for pep in peptides],
            'info':[info_lookup[pep] for pep in peptides]
    })
    elif 'Phenotype' in input_table.columns and 'tissue' not in input_table.columns : 
        # Build a dict to lookup the tissue of the mapped peptides
        #---------------------------------------------------------
        print(f"INFO: {time.ctime()} :: Building a look up table between peptides and phenotype information")
        phenotype_lookup={peptide:phenotype for peptide,phenotype in tqdm(zip(input_table.peptide,input_table.phenotype))}
        print(f"INFO: {time.ctime()} :: Building a look up table between peptides and peptide information")
        info_lookup={peptide:phenotype for peptide,phenotype in tqdm(zip(input_table.peptide,input_table['info']))}
        print(f"INFO: {time.ctime()} :: Building a look up table between peptides and source tissues")
        tissue_lookup={peptide:tissue for peptide,tissue in tqdm(zip(input_table.peptide,input_table.tissue))}
        # Build a table containing the results
        #----------------
        input_peptide=pd.DataFrame({
            'peptide':peptides,
            'allele':HLA,
            'Phenotype':[input_table.loc[input_table.peptide==pep,'phenotype'].to_list()[0] for pep in peptides],
            'info':[input_table.loc[input_table.peptide==pep,'info'].to_list()[0] for pep in peptides],
            'tissue':[input_table.loc[input_table.peptide==pep,'tissue'].to_list()[0] for pep in peptides]
        })
    else: 
         input_peptide=pd.DataFrame({
            'peptide':peptides,
            'allele':HLA})

    ## return the results
    return input_peptide,encoded_input,unmapped
     
def encode_input_for_pia_m(path2input:str, tissue=str)->Tuple[pd.DataFrame,Tuple[np.ndarray,...],pd.DataFrame]:
    """Prepare and input models to PIA-M architecture

    Args:
        path2input (str): The path to load the input table 

    Raises:
        ValueError: Incase the input table is not TAB seperated or has an incorrect number of tables 
        IOError: Incase reading or parsing the input table failed

    Returns:
        Tuple[pd.DataFrame,Tuple[np.ndarray,...],pd.DataFrame]: A tuple of three elements, with the following elements:
        1. pd.DataFrame: a table contain the mapped data  
        2. Tuple[np.ndarray....] a tuple of six tensors representing the encoded arrays (see OmLiT function for more details).
        3. pd.DataFrame: a table contain the unmapped dataset 
    """
     ### Load the input table of peptides and proteins     
     #------------------------------------------------
    try:
        input_table=pd.read_csv(path2input,sep='\t')
        if 'info' in input_table.columns and input_table.shape[1]==3:
            input_table.columns=['peptide','info','allele']
            input_table=input_table[['peptide','allele']]
            input_table['tissue']=[tissue]*input_table.shape[0]
        if input_table.shape[1]==3:
            input_table.columns=['peptide','allele','tissue']
            input_table['tissue'].fillna(tissue, inplace=True)
        elif input_table.shape[1]==2:
            input_table.columns=['peptide','allele']
            input_table['tissue']=[tissue]*input_table.shape[0]
        else:
            raise ValueError(f"In correct number of input columns, excepted 2 or 3 columns, however, the input has: {input_table.shape[0]}")
    except Exception as exp:
        raise IOError(f"{time.ctime()} CRITICAL:: Reading the input table failed with the following error: {exp}")
    ## prepare the inputs
    #--------------------
    ## input to the annotator engine
    input2_OmLiT_annotator=(
        input_table.peptide.to_list(),
        input_table.allele.to_list(),
        input_table.tissue.to_list()
        )
    ## get the encoded and the un-encoded peptides
    encoded_tensors,unmapped_data=linker.annotate_and_encode_input_sequences_no_label(
        input=input2_OmLiT_annotator,max_len=21,proteome=load_proteome(),
        path2cashed_db='/work_ifs/sukmb418/PIA-inference/assets/cashed_database.db',
        path2pseudo_seq='/work_ifs/sukmb418/PIA-inference/assets/pseudosequence.2016.all.X.dat',
        only_one_parent_per_peptide=True 
    )
    unmapped=input_table.loc[unmapped_data[-1]].copy(deep=True).reset_index(drop=True)
    mapped=input_table.drop(unmapped_data[-1],axis=0).reset_index(drop=True)     
    ## return the results
    #--------------------
    return (mapped,
        (np.concatenate([encoded_tensors[0],encoded_tensors[1]],axis=1), 
        np.log10(encoded_tensors[2]+1),
        encoded_tensors[3],
        np.log10(encoded_tensors[4]+1),
        np.log10(encoded_tensors[5]+1)),
        unmapped
    )

def load_proteome()->Dict[str,str]:
     """Load the input FASTA

     Args:
         path2fasta (str): The path to load the FASTA files 

     Returns:
         Dict[str,str]: a dict containing a map between sequence id and sequence names. 
     """
     results=dict()
     for seq in SeqIO.parse('/work_ifs/sukmb418/PIA-inference/assets/filtered_human_sequences_database.fasta','fasta'):
          results[seq.id]=str(seq.seq)
     return results

def map_input_data(sequence, gene_expression, subcellular_location, context, d2g )->Tuple[np.ndarray,...]:
    """Load the input dataset

    Args:
        input_data (data): input tensor obtained by a tensorflow dataset

    Returns:
        Tuple[np.ndarray,...]: A tuple of encoded tensors  
    """
    return {"input_1":sequence,
    "input_2":gene_expression,"input_3":subcellular_location,"input_4":context,"input_5":d2g}

## Define a callback to log the progress
#---------------------------------------
class CustomCallback(tf.keras.callbacks.Callback):
    def on_predict_begin(self, logs=None):
        """will be called at the beging of the predictions. log progress to the .stat file 
        """
        global mapped_peptides
        global log_status_update
        self.num_batches=(mapped_peptides.shape[0]//25_000)+1
        self.log_status_update=log_status_update
        ## print an update message and update the stat file 
        print(f'{time.ctime()} Starting to make predictions using batches, current number of batch is {self.num_batches} \n')
        with open(log_status_update,'w') as writer_stream:
            writer_stream.write(f'0,Starting to make predictions using batches, current number of batch is {self.num_batches} \n')

    def on_predict_batch_begin(self, batch, logs=None):
        """will be called at the beginning of each batch update the stat file with the model progress 

        Args:
            batch (int): the batch index, this input is provided by the fit function. 
            logs (Dict, optional): an optional argument provided by logs to fit the model. Defaults to None.
        """
        print(f'{time.ctime()} Starting to make predictions using batches,running predictions on batch number: {batch} \n')
        with open(self.log_status_update,'w') as writer_stream:
            writer_stream.write(f'{batch/self.num_batches},Running predictions on batch number: {batch} \n')

## Start the main part of the program 
if __name__=='__main__':
    # parse user arguments
    #---------------------
    user_input=parse_user_inputs()
    # define the path
    #----------------
    log_status_update='stat/run.stat'
    ## Log the updates
    print(f'{time.ctime()} loading the model ... \n')
    with open(log_status_update,'w') as writer_stream:
        writer_stream.write(f'0,Loading the model \n')
    model=load_model(user_input['model_index'])
    
    if user_input['model_index']==1:
        # log the update stat
        print(f'{time.ctime()} Encoding and annotating the input ... \n')
        with open(log_status_update,'w') as writer_stream:
            writer_stream.write(f'0,Encoding and annotating the input\n')
        
        ## generate the encoding for PIA-S
        mapped_peptides, encoded_input, unmapped_peptide=encode_input_for_pia_s(user_input['input'])
        
        ## create a dataset from the encoded tensor 
        encoded_input = tf.data.Dataset.from_tensor_slices(encoded_input)
        encoded_input = encoded_input.batch(25_000)

        ## disable AutoShard.
        options = tf.data.Options()
        options.experimental_distribute.auto_shard_policy = tf.data.experimental.AutoShardPolicy.OFF
        encoded_input = encoded_input.with_options(options)
        
        # run prediction 
        predictions=model.predict(encoded_input,verbose=1,callbacks=[CustomCallback()])
        mapped_peptides['predictions']=predictions

    elif user_input['model_index']==2:
        # log the update stat
        print(f'{time.ctime()} Encoding and annotating the input ... \n')
        with open(log_status_update,'w') as writer_stream:
            writer_stream.write(f'0,Encoding and annotating the input\n')
        
        ## load the list of tissues
        #--------------------------
        list_of_tissues=load_list_of_tissues()
        tissue_name=get_tissue_name(user_input['tissue'],list_of_tissues)   

        ## generate the encoding for PIA-M
        mapped_peptides, encoded_tensors, unmapped_peptide=encode_input_for_pia_m(user_input['input'],tissue_name)
        
        ## disable AutoShard.
        options = tf.data.Options()
        options.experimental_distribute.auto_shard_policy = tf.data.experimental.AutoShardPolicy.OFF
        ## build TensorFlow datasets from input tensors
        #----------------------------------------------
        input_seq=tf.data.Dataset.from_tensor_slices(encoded_tensors[0])
        gene_expression=tf.data.Dataset.from_tensor_slices(np.log10(encoded_tensors[1]+1))
        subcellular_location=tf.data.Dataset.from_tensor_slices(encoded_tensors[2])
        context_vector=tf.data.Dataset.from_tensor_slices(np.log10(encoded_tensors[3]+1))
        d2g=tf.data.Dataset.from_tensor_slices(np.log10(encoded_tensors[4]+1))
    
        ## combine all input arrays by zipping the input
        #-----------------------------------------------
        data_feeder=tf.data.Dataset.zip(
                    tuple((input_seq,gene_expression,subcellular_location,context_vector,d2g))
                    ).batch(8192).with_options(options).map(map_input_data)
        predictions=model.predict(data_feeder,verbose=1,callbacks=[CustomCallback()])
        mapped_peptides['predictions']=predictions

    ## write the mapped results, i.e. predictions
    print(f'{time.ctime()} writing the results ... \n')
    with open(log_status_update,'w') as writer_stream:
            writer_stream.write(f'1,Writing the results\n')
    mapped_peptides.to_csv(user_input['output_path'],sep='\t',index=False)
    
    ## write the unmapped results, unmapped input peptides
    print(f'{time.ctime()} writing the unmapped results ... \n')
    with open(log_status_update,'w') as writer_stream:
            writer_stream.write(f'1,Writing the unmapped results\n')
    unmapped_peptide.to_csv(user_input['unmapped_results'],sep='\t',index=False)
    ## update the stat to finished 
    with open(log_status_update,'w') as writer_stream:
            writer_stream.write(f'1,Finished\n')
    ## return 
    sys.exit(0)
