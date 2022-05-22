#!/usr/bin/env 
"""
@author: Hesham ElAbd
@brief: Develop a library of ImFormer models 
"""
## import modules 
#---------------
from typing import *
import tensorflow as tf

## define the ImFormer models
#----------------------------
# implement the transformer layers and architectures 
####################################################
## Note: This implementation is based on 
# https://keras.io/examples/nlp/text_classification_with_transformer/
#########################################

class TransformerBlock(tf.keras.layers.Layer):
    """Creating a transformer layer for sequence classification
    """
    def __init__(self,embed_dim:int,num_heads:int,ff_dim:int,rate=0.1):
        """Create a new Transformer layer 

        Args:
            embed_dim (int): The dimensionality of the embedding layer
            num_heads (int): The number of attention heads 
            ff_dim (int): The number of neurons in the feedforward neural network 
            rate (float, optional): The rate of neural dropout. Defaults to 0.1.

        Returns:
            TransformerBlock: a new transformer block 
        """
        # initialize the base 
        super(TransformerBlock,self).__init__()
        ## add the layer parameters to self so they can be used later with get_cofing and with loading the models 
        self.embed_dim=embed_dim,
        self.num_heads=num_heads
        self.ff_dim=ff_dim
        ## define the compoant of the layer
        self.att = tf.keras.layers.MultiHeadAttention(num_heads=num_heads,key_dim=embed_dim)
        self. ffn = tf.keras.Sequential(
            [tf.keras.layers.Dense(ff_dim,activation="relu"),
            tf.keras.layers.Dense(embed_dim),]
        )
        self.layernorm1=tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2=tf.keras.layers.LayerNormalization(epsilon=1e-6)

        self.dropout1= tf.keras.layers.Dropout(rate)
        self.dropout2=tf.keras.layers.Dropout(rate)
        return 
    
    def get_config(self)->Dict[str,Any]:
        """a helper function to return the parameters of the layer 

        Returns:
            Dict[]
        """

        config=super().get_config()
        config.update({
            'embed_dim':self.embed_dim,
            'num_heads':self.num_heads,
            'ff_dim':self.ff_dim
        })
        return config

    def call(self, inputs, training):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)

class TokenAndPositionEmbedding(tf.keras.layers.Layer):
    """Create an embedding and a position embedding layers"""
    def __init__(self,maxlen:int,vocab_size:int,embed_dim:int):
        """ initialize the position and token embedding 

        Args:
            -maxlen (int): The maximum length of the input sequence 
            -vocab_size (int): The size of the vocabulary 
            -embed_dim (int): The number of dimensions in the embedding layer 

        Returns:
            TokenAndPositionEmbedding: Return a new instance 
        """
        ## initialize the position and token the embeding layer 
        super(TokenAndPositionEmbedding,self).__init__()
        ## saving the parameters to be accessed later by config
        #------------------------------------------------------
        self.maxlen=maxlen
        self.vocab_size=vocab_size
        self.embed_dim=embed_dim
        
        ## now we initialize to embedding layers: 1- Token embedding where we learn a basic embedding, to captsure the meaning of each input fragment
        self.token_emb=tf.keras.layers.Embedding(input_dim=vocab_size, output_dim=embed_dim) # The input dimension --> Number of input dimension; the output dimension is the embedding dimension. This is the meaning of each amino acid.
        self.pos_emb = tf.keras.layers.Embedding(input_dim=maxlen,output_dim=embed_dim) # This the positional embedding 
    
    def call(self,x):
        maxlen=tf.shape(x)[-1]
        positions=tf.range(start=0,limit=maxlen,delta=1) # positions of each input sequence
        positions = self.pos_emb(positions) # encode a position 
        x = self.token_emb(x)
        return x+positions
    
    def get_config(self)->Dict[str,Any]:
        """a helper function to return the parameters of the layer 

        Returns:
            Dict[]
        """
        config=super().get_config()
        config.update({
            'maxlen':self.maxlen,
            'vocab_size':self.vocab_size,
            'embed_dim':self.embed_dim
        })
        return config

## define the modules 
def create_basic_imformer(maxlen:int=27, vocab_size:int=23, 
        embedding_dim:int=4, num_heads:int=2, feedforward:int=64,
        num_blocks:int=1)->tf.keras.Model:
    """Create a transformer model givin user parameters, this is a sequence only model 

    Args:
        maxlen (int, optional): The maximum length of the peptides. Defaults to 27.
        vocab_size (int, optional): The vocabulary size, i.e number of amino acids. Defaults to 23.
        embedding_dim (int, optional): The embedding dimension, i.e. number of amino acid dimensions. Defaults to 4.
        num_heads (int, optional): The number of heads for the multiheaded attension mechanism. Defaults to 2.
        feedforward (int, optional): The number of neurons in the attached block
        num_blocks (int, optional): Number of blocks in the models. Defaults to 1.

    Returns:
        tf.keras.Model: a Transformer based prediction model for predicting peptide HLA interaction
    """
    ## process the input sequence 
    #----------------------------
    inputs=tf.keras.layers.Input(shape=(maxlen,))
    embedding_layer=TokenAndPositionEmbedding(maxlen=maxlen,vocab_size=vocab_size,embed_dim=embedding_dim)
    x=embedding_layer(inputs)
    attention_body=tf.keras.models.Sequential([
        TransformerBlock(embedding_dim,num_heads,feedforward) for _ in range(num_blocks)
        ])
    x=attention_body(x)
    x=tf.keras.layers.GlobalAveragePooling1D()(x)
    x=tf.keras.layers.Dense(20,activation="relu")(x)
    x=tf.keras.layers.Dropout(0.4)(x)
    outputs=tf.keras.layers.Dense(1,activation="sigmoid")(x)
    model=tf.keras.Model(inputs=inputs,outputs=outputs)
    model.compile('adam','binary_crossentropy',
    metrics=["acc",tf.keras.metrics.Recall(),tf.keras.metrics.Precision(),tf.keras.metrics.AUC()])
    return model

##----------------------------------------------------------------
def create_transcription_imformer(maxlen:int=27, vocab_size:int=23, 
        embedding_dim:int=4, num_heads:int=2, feedforward:int=64,
        num_blocks:int=1)->tf.keras.Model:
    """Create transcriptional ImFormer

    Args:
        maxlen (int, optional): [description]. Defaults to 27.
        vocab_size (int, optional): [description]. Defaults to 23.
        embedding_dim (int, optional): [description]. Defaults to 4.
        num_heads (int, optional): [description]. Defaults to 2.
        feedforward (int, optional): [description]. Defaults to 64.
        num_blocks (int, optional): [description]. Defaults to 1.

    Returns:
        tf.keras.Model: [description]
    """
    ## process the input sequence 
    #----------------------------
    inputs=tf.keras.layers.Input(shape=(maxlen,))
    embedding_layer=TokenAndPositionEmbedding(maxlen=maxlen,vocab_size=vocab_size,embed_dim=embedding_dim)
    x=embedding_layer(inputs)
    attention_body=tf.keras.models.Sequential([
        TransformerBlock(embedding_dim,num_heads,feedforward) for _ in range(num_blocks)
        ])
    x=attention_body(x)
    x=tf.keras.layers.GlobalAveragePooling1D()(x)
    x=tf.keras.layers.Dense(20,activation="relu")(x)
    x=tf.keras.layers.Dropout(0.4)(x)
    ## process the input expression  
    #-------------------------------
    transcript_input=tf.keras.layers.Input(shape=(1,))
    y=tf.keras.layers.Dense(1,activation='relu')(transcript_input)
    ## Combine the two inputs
    #------------------------
    combined_stream=tf.keras.layers.concatenate([x,y],axis=-1)
    ## Generate the output model
    #---------------------------  
    outputs=tf.keras.layers.Dense(10,activation='relu')(combined_stream)
    outputs=tf.keras.layers.Dropout(0.2)(outputs)
    outputs=tf.keras.layers.Dense(1,activation="sigmoid")(outputs)
    model=tf.keras.Model(inputs=[inputs,transcript_input],outputs=outputs)
    model.compile('adam','binary_crossentropy',
    metrics=["acc",tf.keras.metrics.Recall(),tf.keras.metrics.Precision(),tf.keras.metrics.AUC()])
    return model

def create_transcription_and_context_imformer(maxlen:int=27, vocab_size:int=23, 
        embedding_dim:int=4, num_heads:int=2, feedforward:int=64,
        num_blocks:int=1)->tf.keras.Model:
    """Create transcriptional ImFormer

    Args:
        maxlen (int, optional): [description]. Defaults to 27.
        vocab_size (int, optional): [description]. Defaults to 23.
        embedding_dim (int, optional): [description]. Defaults to 4.
        num_heads (int, optional): [description]. Defaults to 2.
        feedforward (int, optional): [description]. Defaults to 64.
        num_blocks (int, optional): [description]. Defaults to 1.

    Returns:
        tf.keras.Model: [description]
    """
    ## process the input sequence 
    #----------------------------
    inputs=tf.keras.layers.Input(shape=(maxlen,))
    embedding_layer=TokenAndPositionEmbedding(maxlen=maxlen,vocab_size=vocab_size,embed_dim=embedding_dim)
    x=embedding_layer(inputs)
    attention_body=tf.keras.models.Sequential([
        TransformerBlock(embedding_dim,num_heads,feedforward) for _ in range(num_blocks)
        ])
    x=attention_body(x)
    x=tf.keras.layers.GlobalAveragePooling1D()(x)
    x=tf.keras.layers.Dense(20,activation="relu")(x)
    x=tf.keras.layers.Dropout(0.4)(x)
    ## process the input expression  
    #-------------------------------
    transcript_input=tf.keras.layers.Input(shape=(1,))
    y=tf.keras.layers.Dense(1,activation='relu')(transcript_input)
    ## process the context vector
    #----------------------------
    input_context_vector=tf.keras.layers.Input((16_519,))
    alpha=tf.keras.layers.Dense(1000,activation='relu')(input_context_vector)
    alpha=tf.keras.layers.Dropout(0.4)(alpha)
    alpha=tf.keras.layers.Dense(100,activation='relu')(alpha)
    alpha=tf.keras.layers.Dropout(0.2)(alpha)
    alpha=tf.keras.layers.Dense(10,activation='relu')(alpha)
    ## Combine the two inputs
    #------------------------
    combined_stream=tf.keras.layers.concatenate([x,y,alpha],axis=-1)
    ## Generate the output model
    #---------------------------  
    outputs=tf.keras.layers.Dense(10,activation='relu')(combined_stream)
    outputs=tf.keras.layers.Dropout(0.2)(outputs)
    outputs=tf.keras.layers.Dense(1,activation="sigmoid")(outputs)
    model=tf.keras.Model(inputs=[inputs,transcript_input,input_context_vector],outputs=outputs)
    model.compile('adam','binary_crossentropy',
    metrics=["acc",tf.keras.metrics.Recall(),tf.keras.metrics.Precision(),tf.keras.metrics.AUC()])
    return model

def create_subcellular_location_and_transcription_imformer(maxlen:int=27, vocab_size:int=23, 
        embedding_dim:int=4, num_heads:int=2, feedforward:int=64,
        num_blocks:int=1)->tf.keras.Model:
    """[summary]

    Args:
        maxlen (int, optional): [description]. Defaults to 27.
        vocab_size (int, optional): [description]. Defaults to 23.
        embedding_dim (int, optional): [description]. Defaults to 4.
        num_heads (int, optional): [description]. Defaults to 2.
        feedforward (int, optional): [description]. Defaults to 64.
        num_blocks (int, optional): [description]. Defaults to 1.
        num_subcellular_comps (int, optional): [description]. Defaults to 50.

    Returns:
        tf.keras.Model: [description]
    """
    ## process the input sequence 
    #----------------------------
    inputs=tf.keras.layers.Input(shape=(maxlen,))
    embedding_layer=TokenAndPositionEmbedding(maxlen=maxlen,vocab_size=vocab_size,embed_dim=embedding_dim)
    x=embedding_layer(inputs)
    attention_body=tf.keras.models.Sequential([
        TransformerBlock(embedding_dim,num_heads,feedforward) for _ in range(num_blocks)
        ])
    x=attention_body(x)
    x=tf.keras.layers.GlobalAveragePooling1D()(x)
    x=tf.keras.layers.Dense(20,activation="relu")(x)
    x=tf.keras.layers.Dropout(0.4)(x)
    ## process the input expression level  
    #------------------------------------
    transcript_input=tf.keras.layers.Input(shape=(1,))
    y=tf.keras.layers.Dense(1,activation='relu')(transcript_input)
    ## process the input sub-cellular compartment 
    #--------------------------------------------
    input_subcellular_compartment=tf.keras.layers.Input((1049,))
    z=tf.keras.layers.Dense(100,activation='relu')(input_subcellular_compartment)
    z=tf.keras.layers.Dropout(0.4)(z)
    z=tf.keras.layers.Dense(5,activation='relu')(z)
    ## Combine the three input streams
    #---------------------------------
    combined_stream=tf.keras.layers.concatenate([x,y,z],axis=-1)
    ## Generate the output model
    #---------------------------  
    outputs=tf.keras.layers.Dense(10,activation='relu')(combined_stream)
    outputs=tf.keras.layers.Dropout(0.2)(outputs)
    outputs=tf.keras.layers.Dense(1,activation="sigmoid")(outputs)
    model=tf.keras.Model(inputs=[inputs,transcript_input,input_subcellular_compartment],outputs=outputs)
    model.compile('adam','binary_crossentropy',
    metrics=["acc",tf.keras.metrics.Recall(),tf.keras.metrics.Precision(),tf.keras.metrics.AUC()])
    return model

def create_subcellular_location_transcription_and_contexted_imformer(maxlen:int=27, vocab_size:int=23, 
        embedding_dim:int=4, num_heads:int=2, feedforward:int=64,
        num_blocks:int=1)->tf.keras.Model:
    """[summary]

    Args:
        maxlen (int, optional): [description]. Defaults to 27.
        vocab_size (int, optional): [description]. Defaults to 23.
        embedding_dim (int, optional): [description]. Defaults to 4.
        num_heads (int, optional): [description]. Defaults to 2.
        feedforward (int, optional): [description]. Defaults to 64.
        num_blocks (int, optional): [description]. Defaults to 1.
        num_subcellular_comps (int, optional): [description]. Defaults to 50.

    Returns:
        tf.keras.Model: [description]
    """
    ## process the input sequence 
    #----------------------------
    inputs=tf.keras.layers.Input(shape=(maxlen,))
    embedding_layer=TokenAndPositionEmbedding(maxlen=maxlen,vocab_size=vocab_size,embed_dim=embedding_dim)
    x=embedding_layer(inputs)
    attention_body=tf.keras.models.Sequential([
        TransformerBlock(embedding_dim,num_heads,feedforward) for _ in range(num_blocks)
        ])
    x=attention_body(x)
    x=tf.keras.layers.GlobalAveragePooling1D()(x)
    x=tf.keras.layers.Dense(20,activation="relu")(x)
    x=tf.keras.layers.Dropout(0.4)(x)
    ## process the input expression level  
    #------------------------------------
    transcript_input=tf.keras.layers.Input(shape=(1,))
    y=tf.keras.layers.Dense(1,activation='relu')(transcript_input)
    ## process the input sub-cellular compartment 
    #--------------------------------------------
    input_subcellular_compartment=tf.keras.layers.Input((1049,))
    z=tf.keras.layers.Dense(100,activation='relu')(input_subcellular_compartment)
    z=tf.keras.layers.Dropout(0.4)(z)
    z=tf.keras.layers.Dense(5,activation='relu')(z)
    ## process the context vector
    #----------------------------
    input_context_vector=tf.keras.layers.Input((16_519,))
    alpha=tf.keras.layers.Dense(1000,activation='relu')(input_context_vector)
    alpha=tf.keras.layers.Dropout(0.4)(alpha)
    alpha=tf.keras.layers.Dense(100,activation='relu')(alpha)
    alpha=tf.keras.layers.Dropout(0.4)(alpha)
    alpha=tf.keras.layers.Dense(10,activation='relu')(alpha)
    ## Combine the three input streams
    #---------------------------------
    combined_stream=tf.keras.layers.concatenate([x,y,z,alpha],axis=-1)
    ## Generate the output model
    #---------------------------  
    outputs=tf.keras.layers.Dense(10,activation='relu')(combined_stream)
    outputs=tf.keras.layers.Dropout(0.2)(outputs)
    outputs=tf.keras.layers.Dense(1,activation="sigmoid")(outputs)
    model=tf.keras.Model(inputs=[inputs,transcript_input,input_subcellular_compartment,input_context_vector],
        outputs=outputs)
    model.compile('adam','binary_crossentropy',
    metrics=["acc",tf.keras.metrics.Recall(),tf.keras.metrics.Precision(),tf.keras.metrics.AUC()])
    return model

def create_subcellular_location_transcription_contexted_dist_to_gly_imformer(maxlen:int=27, vocab_size:int=23, 
        embedding_dim:int=4, num_heads:int=2, feedforward:int=64,
        num_blocks:int=1)->tf.keras.Model:
    """[summary]

    Args:
        maxlen (int, optional): [description]. Defaults to 27.
        vocab_size (int, optional): [description]. Defaults to 23.
        embedding_dim (int, optional): [description]. Defaults to 4.
        num_heads (int, optional): [description]. Defaults to 2.
        feedforward (int, optional): [description]. Defaults to 64.
        num_blocks (int, optional): [description]. Defaults to 1.
        num_subcellular_comps (int, optional): [description]. Defaults to 50.

    Returns:
        tf.keras.Model: [description]
    """
    ## process the input sequence 
    #----------------------------
    inputs=tf.keras.layers.Input(shape=(maxlen,))
    embedding_layer=TokenAndPositionEmbedding(maxlen=maxlen,vocab_size=vocab_size,embed_dim=embedding_dim)
    x=embedding_layer(inputs)
    attention_body=tf.keras.models.Sequential([
        TransformerBlock(embedding_dim,num_heads,feedforward) for _ in range(num_blocks)
        ])
    x=attention_body(x)
    x=tf.keras.layers.GlobalAveragePooling1D()(x)
    x=tf.keras.layers.Dense(20,activation="relu")(x)
    x=tf.keras.layers.Dropout(0.4)(x)
    ## process the input expression level  
    #------------------------------------
    transcript_input=tf.keras.layers.Input(shape=(1,))
    y=tf.keras.layers.Dense(1,activation='relu')(transcript_input)
    ## process the input sub-cellular compartment 
    #--------------------------------------------
    input_subcellular_compartment=tf.keras.layers.Input((1049,))
    z=tf.keras.layers.Dense(100,activation='relu')(input_subcellular_compartment)
    z=tf.keras.layers.Dropout(0.4)(z)
    z=tf.keras.layers.Dense(5,activation='relu')(z)
    ## process the context vector
    #----------------------------
    input_context_vector=tf.keras.layers.Input((16_519,))
    alpha=tf.keras.layers.Dense(1000,activation='relu')(input_context_vector)
    alpha=tf.keras.layers.Dropout(0.4)(alpha)
    alpha=tf.keras.layers.Dense(100,activation='relu')(alpha)
    alpha=tf.keras.layers.Dropout(0.2)(alpha)
    alpha=tf.keras.layers.Dense(10,activation='relu')(alpha)
    ## process the distance to glycosylation
    #---------------------------------------
    input_glycosylation_layer=tf.keras.layers.Input((1,))
    beta=tf.keras.layers.Dense(1,activation='relu')(input_glycosylation_layer)
    ## Combine the input streams
    #---------------------------
    combined_stream=tf.keras.layers.concatenate([x,y,z,alpha,beta],axis=-1)
    ## Generate the output model
    #---------------------------  
    outputs=tf.keras.layers.Dense(10,activation='relu')(combined_stream)
    outputs=tf.keras.layers.Dropout(0.2)(outputs)
    outputs=tf.keras.layers.Dense(1,activation="sigmoid")(outputs)
    model=tf.keras.Model(inputs=[inputs,transcript_input,input_subcellular_compartment,
    input_context_vector,input_glycosylation_layer],outputs=outputs)
    model.compile('adam','binary_crossentropy',
    metrics=["acc",tf.keras.metrics.Recall(),tf.keras.metrics.Precision(),tf.keras.metrics.AUC()])
    return model