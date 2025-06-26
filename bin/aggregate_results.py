#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import warnings
#import adata_functions
#from adata_functions import *
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
import ast
import sys
import matplotlib.lines as mlines
from utils import *

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--pipeline_results', type=str, nargs = "+", 
                        help="files containing f1 results with params")                                            
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
 
def main():
    # Parse command line arguments
    args = parse_arguments()
    # Set organism and census_version from arguments
    pipeline_results = args.pipeline_results

    
    predicted_meta_df = pd.DataFrame() 
     
    #for file in os.listdir(pipeline_results):
    for filepath in pipeline_results:
    #filepath = os.path.join(pipeline_results, file)
    # method = filepath.split("/")[-3]
        temp_df = pd.read_csv(filepath,sep="\t")
        # temp_df["method"] = method
        predicted_meta_df = pd.concat([temp_df, predicted_meta_df], ignore_index=True)
     
    organism = predicted_meta_df["organism"].unique()[0]
    # replace "nan" with None
    predicted_meta_df = predicted_meta_df.replace("nan", None)
    #----------weighted f1 results----------------
    # miscellaneous data wrangling
      
    predicted_meta_df["region_match"] = predicted_meta_df.apply(lambda row: row['region'] in row['ref_region'], axis=1)
    predicted_meta_df["reference_acronym"] = predicted_meta_df["reference"].apply(make_acronym)
    predicted_meta_df["reference"] = predicted_meta_df["reference"].str.replace("_", " ")
    predicted_meta_df["study"] = predicted_meta_df["query"].apply(lambda x: x.split("_")[0])
    predicted_meta_df["query"] = predicted_meta_df["query"].str.replace("_", " ")
    
           

    
    predicted_meta_df["disease_state"] = np.where(predicted_meta_df["disease"] == "Control", "Control", "Disease")
    
    
    if organism == "homo_sapiens":
        # data wrangling for missing disease (all controls)
        predicted_meta_df["disease"] = np.where(predicted_meta_df["study"]=="GSE211870", "Control", predicted_meta_df["disease"]) 
    
        # deal with annotation mismatch between gemma queries and curated queries
        predicted_meta_df["dev_stage"] = predicted_meta_df["dev_stage"].apply(map_development_stage) 
        
    # Data wrangling for Rosmap error (dev stage mistakely mapped as "infant")
        predicted_meta_df["dev_stage"] = np.where(predicted_meta_df["study"] == "rosmap" , "late adult", predicted_meta_df["dev_stage"])
        
    # data wrangling for missing Pineda dev stage   
        predicted_meta_df["dev_stage"] = np.where(predicted_meta_df["study"] == "pineda" , "late adult", predicted_meta_df["dev_stage"])

    # data wrangling for Lim sample missing from original metadata
        predicted_meta_df["sex"] = np.where(predicted_meta_df["query"]=="lim_C5382Cd", "male", predicted_meta_df["sex"])
        predicted_meta_df["dev_stage"] = np.where(predicted_meta_df["query"] == "lim_C5382Cd" , "late adult", predicted_meta_df["dev_stage"])


    # data wrangling for sex (Gemmma data uses male:female, conform to this naming scheme)
        predicted_meta_df["sex"] = predicted_meta_df["sex"].str.replace("M", "male")
        predicted_meta_df["sex"] = predicted_meta_df["sex"].str.replace("F", "female")
        # don't know why this is in the data
        predicted_meta_df["sex"] = predicted_meta_df["sex"].str.replace("feM","female")
         
    if organism == "mus_musculus":
        predicted_meta_df["disease_state"] = np.where(predicted_meta_df["disease"].isnull(), "Control", "Disease")
        predicted_meta_df["treatment_state"] = np.where(predicted_meta_df["treatment"].isnull(), "No treatment", "treatment")
        predicted_meta_df["genotype"] = np.where(predicted_meta_df["genotype"].isnull(), "wild type genotype", predicted_meta_df["genotype"])


    predicted_meta_df.to_csv("predicted_meta_combined.tsv", sep="\t", index=False)
        
    # save df for multiqc report
    # get the unique values for each column
    subclass_assignments = (
                            predicted_meta_df
                            .groupby(["subclass", "predicted_subclass"])
                            .size()
                            .unstack(fill_value=0)
                            .reset_index()
                            )

    subclass_assignments.to_csv(os.path.join("subclass_assignments_after_map_mqc.tsv"), sep="\t", index=False)
    
    #do this for class and family
    
    class_assignments = (  predicted_meta_df
                        .groupby(["class", "predicted_class"])  
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    class_assignments.to_csv(os.path.join("class_assignments_after_map_mqc.tsv"), sep="\t", index=False)
    
    family_assignments = ( predicted_meta_df
                        .groupby(["family", "predicted_family"])
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    family_assignments.to_csv(os.path.join("family_assignments_after_map_mqc.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()
    
