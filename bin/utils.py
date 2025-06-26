import pandas as pd
import anndata as ad
import re
import warnings
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
from types import SimpleNamespace


def map_valid_labels(query, ref_keys, mapping_df):
  # deal with differing levels of granularity
  for key in ref_keys:
      print(key)
      original=query[key].unique()
      print(original)
      for og in original:
          # get the highest level in the hierarchy 
          matching_cols = mapping_df.columns[mapping_df.apply(lambda col: og in col.values, axis=0)]
          print(f"Matching columns for {og}: {matching_cols}")
          if len(matching_cols) == 0:
              continue  # likely "unknown", skip
          else:
              level = matching_cols[-1]
              # Check if level is above key in the hierarchy
              if mapping_df.columns.get_loc(level) > mapping_df.columns.get_loc(key):
                  print(f"Level {level} is above level {key} in the hierarchy.")        
                  og_index = query.index[query[key] == og]
                  # Replace the value in "predicted_" column with corresponding predicted value at `level`
                  for idx in og_index:
                      # Find the replacement value from `mapping_df` for this level
                      replacement = query.loc[idx, "predicted_" + level]
                      print(f"Replacing predictions for {og} with {replacement} to match {level}")
                      # replace predicted id with appropriate level
                      query["predicted_" + key] = query["predicted_" + key].astype("object")
                      query.loc[idx, "predicted_" + key] = replacement#.iloc[0]
                      query["predicted_" + key] = query["predicted_" + key].astype("category")

  return query            

def is_correct(adata, ref_keys, mapping_df, level="subclass"):
  # change to string type
  adata.obs["correct_"+level] = adata.obs["predicted_"+level].astype(str) == adata.obs[level].astype(str)
  return adata
  
def stacked_bar_plot(predicted_meta, level="subclass"):
  subclass_assignments = predicted_meta.groupby([level, "correct"]).size().reset_index(name='count')
  pivot_df = subclass_assignments.pivot(index=level, columns='correct', values='count').fillna(0)
  # Sort index (optional)
  pivot_df = pivot_df.sort_index()
  # Plot
  pivot_df.plot(
      kind='barh',
      stacked=True,
      figsize=(20, 20),
      colormap='tab20'  # You can change the color map
  )
  plt.title("Stacked Bar Plot of Correctness")
  plt.xlabel("True Subclass", fontsize=25)
  plt.ylabel("Count", fontsize=25)
  plt.legend(title="Predicted Subclass", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=20, title_fontsize=20)
  plt.tight_layout()
  plt.yticks(fontsize=25)
  plt.xticks(fontsize=25)
  plt.savefig("mapped_correct_boxplot.png", bbox_inches="tight")
  
   
def make_acronym(name):
    # Split on "_" and replace with spaces
    words = name.split("_")
    # Create acronym from the first letter of each word
    acronym = "".join(word[0].upper() for word in words if word)
    return acronym

def map_development_stage(stage):
    # re write dict
    dev_stage_mapping_dict = {
        "HsapDv_0000083": "infant",
        "HsapDv_0000084": "toddler",
        "HsapDv_0000085": "child",
        "HsapDv_0000086": "adolescent",
        "HsapDv_0000088": "adult",
        "HsapDv_0000091": "late adult",
        np.nan: None
    }
    return dev_stage_mapping_dict[stage]
    
def write_factor_summary(df, factors): 
    # summarize the number of unique levels for each item in factors
    # make a value_counts table for each factor

    factor_summary = df[factors].nunique().reset_index()
    factor_summary.columns = ["factor", "levels"]
    factor_summary.to_csv("factor_summary.tsv", sep="\t", index=False) 
    value_counts = df[["study","query","query_region"]].value_counts().reset_index()
    value_counts.to_csv(f"region_study_value_counts.tsv", sep="\t", index=False)