    
#!/user/bin/python3

from pathlib import Path
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
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
import yaml
from utils import *


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--pipeline_results', type=str, nargs = "+", help="Directories containing pipeline results", default=["/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mus_musculus/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/scvi"])                                              
    parser.add_argument('--params_file', type=str, help="Path to the params file", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mus_musculus/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/params.yaml")
    parser.add_argument('--run_name', type=str, help="Name of the original file", default = "ref_50_query_null_cutoff_0_refsplit_dataset_id")
    parser.add_argument('--ref_obs', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mus_musculus/sample/SCT/ref_50_query_null_cutoff_0_refsplit_dataset_id/refs")
    parser.add_argument('--subsample', type=int, default=100, help="Number of cells to subsample from each query dataset")
    parser.add_argument('--ref_keys', type=str, nargs="+", help="Reference keys to map", default=["subclass", "class","family","global"])
    parser.add_argument('--mapping_file', type=str, help="Path to the cell type hierarchy file", default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
    # deal with jupyter kernel arguments
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
def get_ref_regions(ref_obs):
    ref_regions = {}
    for root, dirs, files in os.walk(ref_obs):
        for file in files:
            if file.endswith(".obs.tsv"):
                temp = pd.read_csv(os.path.join(root, file), sep='\t')
                ref_name = file.split(".obs.tsv")[0]
                temp["reference"] = ref_name
                regions = temp["tissue"].unique()
                region="multiple regions" if len(regions) > 1 else regions[0]
                if ref_name not in ref_regions:
                    ref_regions[ref_name] = region
                else:
                    continue
    return ref_regions

def main():
    args = parse_arguments()
    params_file = args.params_file
    pipeline_results = args.pipeline_results
    predicted_meta_df = pd.DataFrame()
    run_name = args.run_name 
    subsample = args.subsample
    ref_obs=args.ref_obs
    ref_keys = args.ref_keys
    mapping_file = args.mapping_file
    reg_regions = get_ref_regions(ref_obs)
    
    mapping_df = pd.read_csv(mapping_file, sep="\t")  # Load the mapping file   
    # Process params_file
    # check if file exists
    if not os.path.exists(params_file):
        raise FileNotFoundError(f"File '{params_file}' not found.")

    with open(params_file, "r") as file:
        parameters_dict = yaml.safe_load(file)  # Parse the YAML file into a Python dictionary

    keys_to_drop = ["ref_collections", "ref_keys", "outdir", 
                    "batch_keys", "relabel_r", "relabel_q", "tree_file","queries_adata"]

    # Use a loop to remove keys
    for key in keys_to_drop:
        parameters_dict.pop(key, None)
    # Optionally, convert the DataFrame back to a dictionary
    predicted_meta_df = pd.DataFrame()       
    # Process tempdf
    for result_path in pipeline_results:  # Assuming tempdf is a list of paths
        method = result_path.split("/")[-1]  # Extract the method from the path
        for root, dirs, files in os.walk(result_path):
            for file in files:
                if "predictions" in file and file.endswith(".tsv"):
                    print(f"Processing file: {file}")
                    study_name = file.split("_")[0] # Extract the study name
                    full_path = os.path.join(root, file)
                    # Extract directory hierarchy
                    ref_name = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname((full_path))))) 
                    query_name = os.path.basename(os.path.dirname(os.path.dirname(full_path)))  # Parent of parent

                    tempdf = pd.read_csv(os.path.join(root, file), sep="\t")  # Read the .tsv file
                    
                    # map valid labels BEFORE subsampling to ensure correct mapping
                    tempdf = map_valid_labels(tempdf, ref_keys, mapping_df)
                    
                    # subsample the DataFrame if it has more than subsample rows
                    if len(tempdf) > subsample:
                        tempdf = tempdf.sample(n=subsample, random_state=42)
                    tempdf["reference"] = ref_name  # Add a reference column
                    tempdf["query"] = query_name
                    tempdf["method"] = method  # Add a method column
                    tempdf["study"] = study_name
                    tempdf["ref_region"] = reg_regions.get(ref_name, "unknown")  # Add reference region
                   
                    # Add the parameters to the DataFrame                 
                    for key, value in parameters_dict.items():
                        tempdf[key] = value
                    predicted_meta_df = pd.concat([predicted_meta_df, tempdf], ignore_index=True)  # Append to the DataFrame
    
    
    predicted_meta_df.to_csv(f"{run_name}_predicted_meta_combined.tsv", sep="\t", index=False)
    
if __name__ == "__main__":
    main()