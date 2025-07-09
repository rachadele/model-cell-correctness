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
from pymer4.models import Lmer

def parse_arguments():
  parser = argparse.ArgumentParser(description="Correct model predictions based on reference keys and mapping file.")
  parser.add_argument('--predicted_meta', type=str, help="Path to the predicted metadata file", default="/space/grp/rschwartz/rschwartz/model_cell_correctness/mus_musculus/results_100/split_by_label/class/Cajal-Retzius_cell/Cajal-Retzius_cell_predicted_meta_subset.tsv")
  parser.add_argument('--mapping_file', type=str, help="Path to the cell type hierarchy file", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
  parser.add_argument('--ref_keys', type=str, nargs="+", help="Reference keys to map", default=["subclass", "class", "family", "global"])
  parser.add_argument('--cell_type', type=str, help="Cell type to analyze", default=None)
  # deal with jupyter kernel arguments
  if __name__ == "__main__":
      known_args, _ = parser.parse_known_args()
      return known_args


def main():
  args = parse_arguments()
  predicted_meta_path = args.predicted_meta
  mapping_file = args.mapping_file
  ref_keys = args.ref_keys
  cell_type = args.cell_type

  # Load predicted metadata
  predicted_meta = pd.read_csv(predicted_meta_path, sep="\t")
  feature_cols = ["outlier_ribo", 
                  "outlier_mito", 
                  "outlier_hb", 
                  "counts_outlier", 
                  "predicted_doublet"]
  formulas = build_formulas_with_doublet_options(predicted_meta, random_effect_column="study")

  for key in ref_keys:
    predicted_meta = is_correct(predicted_meta, level=key)
  outcomes = [f"correct_{key}" for key in ref_keys]

  #for feature in feature_cols:
      #if pd.api.types.is_float_dtype(predicted_meta[feature]):
          ## Melt all outcomes into long-form
          #melted = predicted_meta.melt(
              #id_vars=[feature],
              #value_vars=outcomes,
              #var_name="outcome",
              #value_name="correctness"
          #)
          #plot_feature_violin(melted, feature)
       
      #elif pd.api.types.is_bool_dtype(predicted_meta[feature]) or pd.api.types.is_categorical_dtype(predicted_meta[feature]):
        #melted = predicted_meta.melt(
            #id_vars=[feature],
            #value_vars=outcomes,
            #var_name="outcome",
            #value_name="correctness"
        #)
        
        #plot_feature_boxplot(melted, feature)
        
  models={}
  for formula_name, formula in formulas.items():
    models[formula] = {}
    for key in ref_keys:
        true_col = f"{key}"
        outcome = f"correct_{true_col}"
        new_formula = f"{outcome} ~ {formula}"
        models[formula][key] = {}
        try:
          #fit, df = logistic_regression_on_correctness(predicted_meta, outcome, feature_cols, formula=new_formula)
          model = run_llmer(predicted_meta, new_formula, family="binomial")
        except Exception as e:
          print(f"Error fitting model for {key} with formula {new_formula}: {e}")
          continue
        models[formula][key]["fit"] = model
      
  for formula, dict in models.items():
    coef_list =  []
    for key, subdict in dict.items():
      if "fit" not in subdict:
        print(f"Skipping {key} as fit is not available.")
        continue
      model = subdict["fit"]
      model_df = model.coefs.reset_index().rename(columns={"index": "term"})
      ranef_var= model.ranef_var.reset_index().rename(columns={"index": "term"})
      for row in model_df.iterrows():
        #  if coef_name != "const":  # exclude intercept
        row = row[1].to_dict()
        row["key"] = f"{key}"
        coef_list.append(row)
    # append random effect variance
      for row in ranef_var.iterrows():
        row = row[1].to_dict()
        row["Estimate"] = row["Var"]
        row["2.5_ci"] = row["Var"] - 2 * row["Std"]
        row["97.5_ci"] = row["Var"] + 2 * row["Std"]
        row["key"] = f"{key}"
        coef_list.append(row)
    # Combine into a single DataFrame
    coef_df = pd.DataFrame(coef_list)
    
    # save coef df to file
    if cell_type is None:
      cell_type = "all"
    coef_df.to_csv(f"{formula.replace(' ', '_')}_{cell_type}_coefficients.tsv", sep="\t", index=False)
  # Plot coefficients 
    plot_coefficients(coef_df, formula = formula.replace(" ", "_"))

  
if __name__ == "__main__":
  main() 