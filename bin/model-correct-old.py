from utils import *
import numpy as np
import pandas as pd
from pathlib import Path
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import statsmodels.api as sm
import statsmodels.formula.api as smf
from utils import *
from collections import defaultdict

def parse_arguments():
  parser = argparse.ArgumentParser(description="Correct model predictions based on reference keys and mapping file.")
  parser.add_argument('--predicted_meta', type=str, help="Path to the predicted metadata file", default="/space/grp/rschwartz/rschwartz/model_cell_correctness/mus_musculus/results_100/aggregated_results/predicted_meta_combined.tsv")
  parser.add_argument('--mapping_file', type=str, help="Path to the cell type hierarchy file", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
  parser.add_argument('--ref_keys', type=str, nargs="+", help="Reference keys to map", default=["subclass", "class", "family", "global"])

  # deal with jupyter kernel arguments
  if __name__ == "__main__":
      known_args, _ = parser.parse_known_args()
      return known_args


def main():
  args = parse_arguments()
  predicted_meta_path = args.predicted_meta
  mapping_file = args.mapping_file
  ref_keys = args.ref_keys
  
  # Load predicted metadata
  predicted_meta = pd.read_csv(predicted_meta_path, sep="\t")
  
  # subset to only cutoff=0, reference= whole cortex, method= scvi
  predicted_meta = predicted_meta[(predicted_meta["cutoff"] == 0) &
                                  (predicted_meta["reference"] == "whole cortex") &
                                  (predicted_meta["method"] == "scvi")]
 
  random_effect_column = "study" 
  # logistic regression of correct vs incorrect predictions
  feature_cols = ["outlier_ribo", 
                  "outlier_mito", 
                  "outlier_hb", 
                  "counts_outlier", 
                  "doublet_score",
                  "predicted_doublet"]
 # method and reference removed
  
  formulas = ["outlier_ribo + outlier_mito + outlier_hb + counts_outlier + doublet_score",
              "outlier_ribo + outlier_mito + outlier_hb + counts_outlier + predicted_doublet"]

  # make default dict
  
  fits = {}

  for formula in formulas:
    fits[formula] = {}
    for key in ref_keys:
      fits[formula][key] = {}
    # preditcted_meta_subset = predicted_meta[predicted_meta["key"] == key]
      true_col = f"{key}"
      predicted_meta = is_correct(predicted_meta, level=key)
      outcome = f"correct_{true_col}"
      fit, df = logistic_regression_on_correctness(predicted_meta, outcome, feature_cols, formula=formula)
      fits[formula][key]["fit"] = fit
     # fits[formula][key]["predicted_meta"] = df
    
  for formula, dict in fits.items():

    coef_list =  []
    for key, subdict in dict.items():
      print(key) 
      fit = subdict["fit"]
      for coef_name, coef_value in fit.params.items():
        #  if coef_name != "const":  # exclude intercept
        coef_list.append({
            "key": key,
            "feature": coef_name,
            "coefficient": coef_value
        })
    # Combine into a single DataFrame
    coef_df = pd.DataFrame(coef_list)
  # Plot coefficients 
    plot_coefficients(coef_df, out_prefix = formula.replace(" ", "_"))

  # Assume ref_keys and predicted_meta are defined
  outcomes = [f"correct_{key}" for key in ref_keys]


  for feature in feature_cols:
      if pd.api.types.is_float_dtype(predicted_meta[feature]):
          # Melt all outcomes into long-form
          melted = predicted_meta.melt(
              id_vars=[feature],
              value_vars=outcomes,
              var_name="outcome",
              value_name="correctness"
          )

          g = sns.catplot(
              data=melted,
              x="correctness",
              y=feature,
              col="outcome",
              kind="violin",
              inner="quartile",
              height=4,
              aspect=0.8,
              hue="correctness",
              palette="muted",
              col_wrap=4,
              sharey=True
          )
          g.set_titles("{col_name}")
          g.fig.subplots_adjust(top=0.9)
          g.fig.suptitle(f"{feature} by correctness (all outcomes)")
          g.tight_layout()

          outdir = "combined_plots"
          os.makedirs(outdir, exist_ok=True)
          g.savefig(os.path.join(outdir, f"{feature}_by_correctness_all.png"))

      elif pd.api.types.is_bool_dtype(predicted_meta[feature]) or pd.api.types.is_categorical_dtype(predicted_meta[feature]):
        melted = predicted_meta.melt(
            id_vars=[feature],
            value_vars=outcomes,
            var_name="outcome",
            value_name="correctness"
        )

        g = sns.catplot(
            data=melted,
            x=feature,
            hue="correctness",
            col="outcome",
            kind="count",
            col_wrap=4,
            height=4,
            aspect=1.2,
            sharey=False
        )
        g.set_titles("{col_name}")
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle(f"{feature} counts by correctness (all outcomes)")
        g.set_xticklabels(rotation=45)

        outdir = "combined_plots"
        os.makedirs(outdir, exist_ok=True)
        g.savefig(os.path.join(outdir, f"{feature}_by_correctness_all.png"))

  
if __name__ == "__main__":
  main() 