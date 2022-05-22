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
from tqdm import tqdm 
from Bio import SeqIO
## define the analysis function 
#--------------------------
def parse_user_inputs()->Dict[str,str]:
     """Parse use input and returns a dict containing user parameters 

     Returns:
         Dict[str,str]: a dict containing user provided elements 
     """
     # Define the parse to parse the results
     #--------------------------------------
     parser=argparse.ArgumentParser(description="Predicting the interaction between an input list of peptides and a list of HLA proteins using either PIA-S or PIA-M")
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
     default=os.path.join(os.getcwd(),'prediction_results.tsv'),
     type=str
     )
     parser.add_argument('--unmapped_results',
     help="The path to write the unmapped results, defaults to unmapped_results.tsv in the current working directory.",
     default=os.path.join(os.getcwd(),'unmapped_results.tsv'),
     type=str
     )
     ## Parse input parameters
     #-----------------------
     args=parser.parse_args()
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
     if args.model_index not in  [1,2,3,4]:
          sys.stderr.write(f"{time.ctime()} CRITICAL:: The model index is not valid, the following four models are only valid: 1. is PIA-S trained on public databases, 2. PIA-S trained on public and inhouse datasets, 3. is PIA-M trained on public databases and 4.\
     is PIA-M trained on all public and inhouse databases. while you input is: {args.model_index}\n")
     else: 
          model_index=args.model_index
     # Check the validity of the output result file 
     if os.access(os.path.dirname(args.output_path), os.W_OK):
          results_path=args.output_path
     else:
          os.stderr.write(f"{time.ctime()} CRITICAL:: invalid writting path to write the output results, user: {os.getlogin()} does not not have writting access at: {args.output_path}, are you sure the output path is valid.")
     # Check the validity of the output result file 
     if os.access(os.path.dirname(args.unmapped_results), os.W_OK):
          unmapped_results=args.unmapped_results
     else:
          os.stderr.write(f"{time.ctime()} CRITICAL:: invalid writting path to write the unmapped results, user: {os.getlogin()} does not not have writting access at: {args.unmapped_results}, are you sure the output path is valid.")
     # return the results
     #-------------------
     return {
          'input':input_path,
          'model_index':model_index,
          'output_path':results_path,
          'unmapped_results':unmapped_results
     }

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
          if model_index==0:
               model=create_basic_imformer(maxlen=55,vocab_size=28,embedding_dim=32,num_heads=4,feedforward=64,num_blocks=3)
               model.load_weights('~helabd/PIA/assets/PIA_Full_public_only.h5')
               print(f"{time.ctime()} INFO:: PIA has been loaded")
          elif model_index==1:
               model=create_basic_imformer(maxlen=55,vocab_size=28,embedding_dim=32,num_heads=4,feedforward=64,num_blocks=3)
               model.load_weights('~helabd/PIA/assets/PIA_Full_public_and_inhouse.h5')
               print(f"{time.ctime()} INFO:: PIA has been loaded")
          elif model_index==2:
               model=create_subcellular_location_transcription_contexted_dist_to_gly_imformer(maxlen=55,vocab_size=28,embedding_dim=32,num_heads=4,feedforward=64,num_blocks=3) 
               model.load_weights('~helabd/PIA/assets/PIA_M_public.h5')
               print(f"{time.ctime()} INFO:: PIA has been loaded")
          elif model_index==3:
               model=create_subcellular_location_transcription_contexted_dist_to_gly_imformer(maxlen=55,vocab_size=28,embedding_dim=32,num_heads=4,feedforward=64,num_blocks=3) 
               model.load_weights('~helabd/PIA/assets/PIA_M_public_and_inhouse.h5')
               print(f"{time.ctime()} INFO:: PIA has been loaded")
          else: 
               raise ValueError(f"{time.ctime()}:: Unsupported model index, current model only support 4 indices,\
                    1 --> PIA-S (public), 2 --> PIA-S (public+inhouse), 3 --> PIA-M (public), 4 --> PIA-M (public+inhouse). while your input is: {model_index}")
     return model

def encode_input_for_pia_s(path2input:str, model_index:int)->Tuple[pd.DataFrame,np.ndarray,pd.DataFrame]:
     """prepare and encode the data for running inferences using PIA-S

     Args:
         path2input (str): The path to load the input file 
         model_index (int): The index of the model to use

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
          if input_table.shape[1]==3:
               input_table.columns=['peptide','allele','tissue']
          elif input_table.shape[1]==2:
               input_table.columns=['peptide','allele']
          else:
               raise ValueError(f"In correct number of input columns, excepted 2 or 3 columns, however, the input has: {input_table.shape[0]}")
     except Exception as exp:
          raise IOError(f"{time.ctime()} CRITICAL:: Reading the input table failed with the following error: {exp}")
     
     ### Load the pseudo sequence
     #---------------------------
     pseudo_sequences=pd.read_csv('~helabd/PIA/assets/pseudosequence.2016.all.X.dat',header=None,sep='\t')
     pseudo_sequences.columns=['Allele','pseudo_sequence']

     ### get unmapped due-to un matached alleles
     #------------------------------------------
     pseudo_sequences_look_up=pd.read_csv('~helabd/PIA/assets/pseudosequence.2016.all.X.dat',header=None,sep='\t')
     pseudo_sequences_look_up.columns=['Allele','pseudo_sequence']
     ### Remove un filtered 
     print(f"{time.ctime()} INFO:: filtering the input table with {input_table.shape[0]} elements.")
     unmapped=input_table.loc[~((input_table.peptide.str.isalpha()) & (input_table.allele.isin(pseudo_sequences_look_up.Allele))),]
     mapped=input_table.loc[((input_table.peptide.str.isalpha()) & (input_table.allele.isin(pseudo_sequences_look_up.Allele))),]
     print(f"{time.ctime()} INFO:: finished with filtering the inputs, number of unmapped is: {unmapped.shape[0]} while number of mapped is: {mapped.shape[0]}")
     
     ### Encoding input peptides:
     #---------------------------
     print(f"{time.ctime()} INFO:: Extracted the pseudo-sequences for each allele")
     input_to_PIA_S, peptides, HLA=[],[],[]
     for allele_name, grouped_peptides in tqdm(mapped.groupby('allele')):
          pseudo_seq=pseudo_sequences_look_up.loc[pseudo_sequences_look_up.Allele==allele_name,'Allele'].to_list()[0]
          input_to_PIA_S.extend([peptide+pseudo_seq for peptide in grouped_peptides.peptide])
          peptides.extend(grouped_peptides.peptide)
          HLA.extend(grouped_peptides.allele)
     ## Load the Tokenizer
     #--------------------
     if model_index==1:
          with open('~helabd/PIA/assets/tokenizer_PIA_S_public.pickle','rb') as buffer_reader:
               tokenizer=pickle.load(buffer_reader)
     elif model_index==2:
          with open('~helabd/PIA/assets/tokenizer_PIA_S_public_and_inhouse.pickle','rb') as buffer_reader:
               tokenizer=pickle.load(buffer_reader)
     ### Encode the input
     #-------------------
     encoded_input=tf.keras.preprocessing.sequence.pad_sequences(
            sequences=tokenizer.texts_to_sequences(input_to_PIA_S),
            dtype=np.int32,maxlen=55,padding="pre")
     ### Make a dataframe from the results
     #------------------------------------
     input_peptide=pd.DataFrame({
          'peptide':peptides,
          'allele':HLA,
     })
     ## return the results
     return input_peptide,encoded_input,unmapped
     
def encode_input_for_pia_m(path2input:str)->Tuple[pd.DataFrame,Tuple[np.ndarray,...],pd.DataFrame]:
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
          input_table=pd.read_csv(path2input,sep='\t').iloc[:,:2]
          if input_table.shape[1]==3:
               input_table.columns=['peptide','allele','tissue']
               input_table['tissue'].fillna('total PBMC', inplace=True)
          elif input_table.shape[1]==2:
               input_table.columns=['peptide','allele']
               input_table['tissue']=['total PBMC']*input_table.shape[0]
          else:
               raise ValueError(f"In correct number of input columns, excepted 2 or 3 columns, however, the input has: {input_table.shape[0]}")
     except Exception as exp:
          raise IOError(f"{time.ctime()} CRITICAL:: Reading the input table failed with the following error: {exp}")
     ## prepare the inputs
     #--------------------
     ## input to the annotator engine
     input2_OmLiT_annotator=(
          input_table.peptides.to_list(),
          input_table.allele.to_list(),
          input_table.tissue.to_list()
          )
     ## get the encoded and the un-encoded peptides
     encoded_tensors,unmapped_data=linker.annotate_and_encode_input_sequences_no_label(
          input=input2_OmLiT_annotator,max_len=21,proteome=load_proteome(),
          path2cashed_db='~helabd/PIA/assets/cashed_database.db',
          path2pseudo_seq='~helabd/PIA/assets/pseudosequence.2016.all.X.dat',
          only_one_parent_per_peptide=True 
     )
     unmapped=input_table.loc[unmapped_data[-1]].copy(deep=True).reset_index(drop=True)
     mapped=input_table.drop(unmapped_data[-1],axis=0).reset_index(drop=True)     
     ## return the results
     #--------------------
     return (mapped,
               (np.concatenate([encoded_tensors[0],encoded_tensors[1]],axis=1), 
               np.log10(encoded_tensors[2]+1),encoded_tensors[3], np.log10(encoded_tensors[4]+1),np.log10(encoded_tensors[5]+1)),
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
     for seq in SeqIO.parse('~helabd/PIA/assets/filtered_human_sequences_database.fasta','fasta'):
          results[seq.id]=str(seq.seq)
     return results

## Start the main part of the program 
if __name__=='__main__':
     user_input=parse_user_inputs()
     model=load_model(user_input['model_index'])
     if user_input['model_index']==1 or user_input['model_index']==2:
          mapped_peptides, encoded_input, unmapped_peptide=encode_input_for_pia_s(user_input['input'],user_input['model_index'])
          predictions=model.predict(encoded_input,batch_size=10_000,verbose=1)
          mapped_peptides['predictions']=predictions
     elif user_input['model_index']==3 or user_input['model_index']==4:
          mapped_peptides, encoded_input, unmapped_peptide=encode_input_for_pia_m(user_input['input'])
          predictions=model.predict(encoded_input,batch_size=10_000,verbose=1)
          mapped_peptides['predictions']=predictions
     ### Write the mapped results
     mapped_peptides.to_csv(user_input['output_path'],sep='\t',index=False)
     unmapped_peptide.to_csv(user_input['unmapped_results'],sep='\t',index=False)
