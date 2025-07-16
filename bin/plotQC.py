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
    parser.add_argument('--predicted_meta', type=str, default="/space/grp/rschwartz/rschwartz/model_cell_correctness/mus_musculus/results_10/aggregated_results/predicted_meta_combined.tsv",
                        help="Path to the predicted metadata file")
    parser.add_argument("--ref_keys", type=str, nargs="+",
                        help="Reference keys to map",
                        default=["subclass", "class", "family", "global"]) 


    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
 
def main():
    # Parse command line arguments
    args = parse_arguments()
    ref_keys = args.ref_keys
    predicted_meta = args.predicted_meta
    predicted_meta_df = pd.read_csv(predicted_meta, sep="\t")
       # filter out results with cutoff above cutoff
       
    #predicted_meta_df = predicted_meta_df[predicted_meta_df["cutoff"] <= cutoff]

    #predicted_meta_df.to_csv(f"predicted_meta_combined_cutoff_{cutoff}.tsv", sep="\t", index=False)
       
    multiqc_dir = "multiqc"
    os.makedirs(multiqc_dir, exist_ok=True) 
    # save df for multiqc report
    # get the unique values for each column
    subclass_assignments = (
                            predicted_meta_df
                            .groupby(["subclass", "predicted_subclass"])
                            .size()
                            .unstack(fill_value=0)
                            .reset_index()
                            )

    subclass_assignments.to_csv(os.path.join(multiqc_dir,"pred_subclass_assignments_mqc.tsv"), sep="\t", index=False)
    
    #do this for class and family
    
    class_assignments = (  predicted_meta_df
                        .groupby(["class", "predicted_class"])  
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    class_assignments.to_csv(os.path.join(multiqc_dir, "pred_class_assignments_mqc.tsv"), sep="\t", index=False)
    
    family_assignments = ( predicted_meta_df
                        .groupby(["family", "predicted_family"])
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    family_assignments.to_csv(os.path.join(multiqc_dir, "pred_family_assignments_mqc.tsv"), sep="\t", index=False)
    
    if "global" in ref_keys:
        global_assignments = ( predicted_meta_df
                    .groupby(["global", "predicted_global"])
                    .size()
                    .unstack(fill_value=0)
                    .reset_index()
                    )
        global_assignments.to_csv(os.path.join(multiqc_dir, "pred_global_assignments_mqc.tsv"), sep="\t", index=False)


    # breakdown correctness by treatment, sex, development stage, query, reference
    if "treatment" in predicted_meta_df.columns:
        treatment_correctness = (predicted_meta_df
                            .groupby(["treatment", "correct_subclass"])
                            .size()
                            .unstack(fill_value=0)
                            .reset_index()
                            )
        treatment_correctness.to_csv(os.path.join(multiqc_dir, "treatment_correctness_mqc.tsv"), sep="\t", index=False)
    
    if "genotype" in predicted_meta_df.columns:
        genotype_correctness = (predicted_meta_df
                            .groupby(["genotype", "correct_subclass"])
                            .size()
                            .unstack(fill_value=0)
                            .reset_index()
                            )
        genotype_correctness.to_csv(os.path.join(multiqc_dir, "genotype_correctness_mqc.tsv"), sep="\t", index=False)
    
    sex_correctness = (predicted_meta_df
                        .groupby(["sex", "correct_subclass"])
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    sex_correctness.to_csv(os.path.join(multiqc_dir, "sex_correctness_mqc.tsv"), sep="\t", index=False)
    
    study_correctness = (predicted_meta_df
                        .groupby(["study", "correct_subclass"])
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    study_correctness.to_csv(os.path.join(multiqc_dir, "study_correctness_mqc.tsv"), sep="\t", index=False)
    
    ref_correctness = (predicted_meta_df
                        .groupby(["reference", "correct_subclass"])
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    ref_correctness.to_csv(os.path.join(multiqc_dir, "ref_correctness_mqc.tsv"), sep="\t", index=False)
    
    method_correctness = (predicted_meta_df
                        .groupby(["method", "correct_subclass"])
                        .size()
                        .unstack(fill_value=0)
                        .reset_index()
                        )
    method_correctness.to_csv(os.path.join(multiqc_dir, "method_correctness_mqc.tsv"), sep="\t", index=False)

    #if "age" in predicted_meta_df.columns:
      #  #check if float
       # if predicted_meta_df["age"].dtype == "float64":
            
if __name__ == "__main__":
    main()
    
