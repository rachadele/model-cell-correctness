from utils import *
import numpy as np
import pandas as pd
from pathlib import Path
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import statsmodels.api as sm


def parse_arguments():
  parser = argparse.ArgumentParser(description="Correct model predictions based on reference keys and mapping file.")
  parser.add_argument('--predicted_meta', type=str, help="Path to the predicted metadata file", default="/space/grp/rschwartz/rschwartz/model_cell_correctness/results_10/aggregated_results/predicted_meta_combined.tsv")
  parser.add_argument('--mapping_file', type=str, help="Path to the cell type hierarchy file", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
  parser.add_argument('--ref_keys', type=str, nargs="+", help="Reference keys to map", default=["subclass", "class", "family", "global"])

  # deal with jupyter kernel arguments
  if __name__ == "__main__":
      known_args, _ = parser.parse_known_args()
      return known_args


def logistic_regression_on_correctness(predicted_meta, pred_col, true_col, feature_cols):
  for col in feature_cols:
      if pd.api.types.is_bool_dtype(predicted_meta[col]):
          predicted_meta[col] = predicted_meta[col].astype(int)

  predicted_meta[f"correct_{true_col}"] =   predicted_meta[f"correct_{true_col}"].astype(int) 
  X = sm.add_constant(predicted_meta[feature_cols])
  y = predicted_meta[f"correct_{true_col}"]

  result = sm.Logit(y, X).fit()
  return result, predicted_meta

def plot_coefficients(coef_df):
    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=coef_df,
        y="feature",
        x="coefficient",
        hue="key",
        orient="h"
    )
    plt.axvline(0, color="gray", linestyle="--")
    plt.title("Logistic Regression Coefficients by Key")
    plt.xlabel("Coefficient (Log-Odds)")
    plt.ylabel("Feature")
    plt.tight_layout()
    plt.savefig("logistic_regression_coefficients.png")

def main():
  args = parse_arguments()
  predicted_meta_path = args.predicted_meta
  mapping_file = args.mapping_file
  ref_keys = args.ref_keys
  
  # Load predicted metadata
  predicted_meta = pd.read_csv(predicted_meta_path, sep="\t")
  # Load mapping file
  mapping_df = pd.read_csv(mapping_file, sep="\t")
  
  # logistic regression of correct vs incorrect predictions
  feature_cols = ["outlier_ribo","outlier_mito", "outlier_hb", "counts_outlier", "doublet_score"]
  fits = {key: {"fit": None, "predicted_meta": None} for key in ref_keys}
  for key in ref_keys:
   # preditcted_meta_subset = predicted_meta[predicted_meta["key"] == key]
    pred_col = f"predicted_{key}"
    true_col = f"{key}"
    predicted_meta = is_correct(predicted_meta, level=key)
    fit, df = logistic_regression_on_correctness(predicted_meta, pred_col, true_col, feature_cols)
    fits[key]["fit"] = fit
    fits[key]["predicted_meta"] = df
    # SUMMARIZE FIT
    
  coef_list = []

  for key in ref_keys:
      fit = fits[key]["fit"]
      for coef_name, coef_value in fit.params.items():
          if coef_name != "const":  # exclude intercept
              coef_list.append({
                  "key": key,
                  "feature": coef_name,
                  "coefficient": coef_value
              })

  # Combine into a single DataFrame
  coef_df = pd.DataFrame(coef_list)
  # Plot coefficients 
  plot_coefficients(coef_df)



    
if __name__ == "__main__":
  main()
