#!/user/bin/python3

from pathlib import Path
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import warnings
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import statsmodels as sm
from scipy import stats
import matplotlib.pyplot as plt
set_seed = 42
import random
import numpy as np
random.seed(42)


# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--predicted_meta', type=str, help="predicted metadata", default = "/space/grp/rschwartz/rschwartz/model_cell_correctness/mus_musculus/results_100/aggregated_results/predicted_meta_combined.tsv")   
  parser.add_argument("--ref_keys", type=str, nargs="+", help="Reference keys to map", default=["subclass", "class", "family", "global"])
  # deal with jupyter kernel arguments
  if __name__ == "__main__":
      known_args, _ = parser.parse_known_args()
      return known_args
  
  return parser.parse_args()


def main():
  args = parse_arguments()
  predicted_meta = args.predicted_meta
	
	# Read the label f1 results
  df = pd.read_csv(predicted_meta, sep="\t")
	# split by key and label and save to individual files
 # for key, group in df.groupby("key"):
  #  for label, sub_group in group.groupby("label"):
  levels = args.ref_keys
  for key in levels:
    for label in df[key].unique():
      # filter the dataframe for the current key and label
      label_df = df[(df[key] == label)]
      # check if the label_df is empty
      if label_df.empty:
        print(f"No data for key: {key}, label: {label}")
        continue
      # check if the label_df has more than 1 row
      if len(label_df) < 10:
        print(f"Not enough data for key: {key}, label: {label}")
        continue
  
    # plot the f1 score distribution
      label_name = label.replace(" ", "_").replace("/", "_")
      os.makedirs(f"{key}/{label_name}", exist_ok=True)
      output_file = f"{key}/{label_name}/{label_name}_predicted_meta_subset.tsv"
      label_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
	main()