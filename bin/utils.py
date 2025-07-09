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
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
from pymer4.models import Lmer

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

def is_correct(df, level="subclass"):
  # change to string type
  df["correct_"+level] = df["predicted_"+level].astype(str) == df[level].astype(str)
  return df
  
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
        "nan": None,
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
    
def logistic_regression_on_correctness(predicted_meta, outcome, feature_cols, formula=None):
    df = predicted_meta.copy()
    # Convert boolean to int, and ensure numerics are floats
    for col in feature_cols:
        if pd.api.types.is_bool_dtype(df[col]):
            df[col] = df[col].astype(int)
        elif pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].astype(float)
        # string/categorical will be handled by formula
    
    df[outcome] = df[outcome].astype(int)
    terms = []
    # Fit logistic regression
    fit = smf.logit(formula=formula, data=df).fit()
    return fit, df

def run_llmer(df, formula, family="binomial"):
  print(f"Fitting model with formula: {formula} and family: {family}")
  model = Lmer(formula, data=df, family=family)
  #model.set_factors(["study"])
  model.fit()
  return(model)

def plot_random_slopes(model, cell_type):
 # Get all fixed effect names (from .coefs index)
  all_terms = model.coefs.index.tolist()
  for term in all_terms:
      try:
          print(f"Plotting random slopes for: {term}")
          model.plot(param=term, figsize=(8, 6))
          plt.title(f"Effect of '{term}' by Study")
          plt.tight_layout()
          plt.savefig(f"{cell_type}_random_slope_{term}.png")
          plt.close()
      except Exception as e:
          print(f"Skipping {term}: {e}")


def plot_coefficients(coef_df, formula="logistic_regression"):
    os.makedirs(formula, exist_ok=True)
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=coef_df,
        y="term",
        x="Estimate",
        hue="key",
        orient="h"
    )
    # add confidence intervals
    for bar, (_,row) in zip(ax.patches, coef_df.iterrows()):
        bar_center = bar.get_y() + bar.get_height() / 2
        ci_low = row["2.5_ci"]
        ci_high = row["97.5_ci"]
        estimate = row["Estimate"]
        
        plt.errorbar(
            x=estimate,
            y=bar_center,
            xerr=[[estimate - ci_low], [ci_high - estimate]],
            fmt='o',
            color='black',
            capsize=5,
            elinewidth=1,
            markersize=3
        )
        
        
    plt.axvline(0, color="gray", linestyle="--")
    formula_str = formula.replace("_", " ")
    plt.title(f"{formula_str} Coefficients")
    plt.xlabel("Log Odds")
    plt.ylabel("Feature")
    plt.tight_layout()
    plt.savefig(os.path.join(formula,"logistic_regression_coefficients.png"))
    
def plot_feature_violin(melted, feature):
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
  g.fig.subplots_adjust(top=0.8)
  g.fig.suptitle(f"{feature} by correctness")
  g.tight_layout()

  outdir = "combined_plots"
  os.makedirs(outdir, exist_ok=True)
  g.savefig(os.path.join(outdir, f"{feature}_by_correctness.png"))


def plot_feature_boxplot(melted, feature):
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
  g.fig.subplots_adjust(top=0.8)
  g.fig.suptitle(f"{feature} counts by correctness")
  g.set_xticklabels(rotation=45)

  outdir = "combined_plots"
  os.makedirs(outdir, exist_ok=True)
  g.savefig(os.path.join(outdir, f"{feature}_by_correctness.png"))



def build_formulas_with_doublet_options(predicted_meta, random_effect_column="study"):
  """
  Build two model formulas: one using predicted_doublet, one using doublet_score.
  Only includes features with more than one unique value.
  
  Returns:
      formulas (dict): keys are 'predicted_doublet' and 'doublet_score'
  """
  base_features = ["outlier_ribo", "outlier_mito", "outlier_hb", "counts_outlier"]
  doublet_features = ["predicted_doublet", "doublet_score"]
  
  formulas = {}
  for doublet_col in doublet_features:
    feature_cols = base_features + [doublet_col]
    valid_feature_cols = []

    for col in feature_cols:
      if col in predicted_meta.columns:
        unique_vals = predicted_meta[col].dropna().unique()
        if len(unique_vals) > 1:
          valid_feature_cols.append(col)
        else:
          print(f"Skipping {col}: only one unique value.")
      else:
        print(f"Skipping {col}: column not found.")

    if len(valid_feature_cols) > 0:
      formula = " + ".join(valid_feature_cols) + f" + (1 | {random_effect_column})"
      formulas[doublet_col] = formula
    else:
      print(f"No valid features found for formula with {doublet_col}.")

  if not formulas:
    raise ValueError("No valid formulas could be constructed.")
  return formulas