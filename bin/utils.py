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
  df["correct_" + level] = df["predicted_" + level].astype(str) == df[level].astype(str)
  return df

   
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
    
def logistic_regression_on_correctness(predicted_meta, outcome, formula=None):
    df = predicted_meta.copy()
    # Convert boolean to int, and ensure numerics are floats 
    df[outcome] = df[outcome].astype(int)
    # Fit logistic regression
    fit = smf.logit(formula=formula, data=df).fit()
    return fit

def run_llmer(df, formula, family="binomial"):
  """Fit a binomial GLMM via direct rpy2 → lme4::glmer.

  Returns a SimpleNamespace with `.coefs` (term × {Estimate, SE, Z-stat, P-val,
  2.5_ci, 97.5_ci, OR/Prob and CIs, Sig}) and `.ranef_var` (grp × {Var, Std}),
  matching the attributes expected by model_correct.py. We bypass pymer4 because
  its post-fit summary extraction returns NULL on large dataframes in this env
  (pandas2R quirk); direct glmer fits the same data without issue.
  """
  print(f"Fitting model with formula: {formula} and family: {family}")
  if family != "binomial":
    raise NotImplementedError(f"only binomial family supported, got {family}")
  from rpy2 import robjects
  from rpy2.robjects import pandas2ri, conversion, default_converter
  from rpy2.robjects.packages import importr
  cv = default_converter + pandas2ri.converter
  with conversion.localconverter(cv):
    importr("lme4")
    robjects.globalenv["dat"] = pandas2ri.py2rpy(df)
    robjects.r(f'fit_ <- glmer({formula}, data=dat, family=binomial)')
    coef_mat = np.array(robjects.r('coef(summary(fit_))'))
    rownames = list(robjects.r('rownames(coef(summary(fit_)))'))
    vc = robjects.r('as.data.frame(VarCorr(fit_))')
    if not isinstance(vc, pd.DataFrame):
      vc = pandas2ri.rpy2py(vc)
  coefs = pd.DataFrame(coef_mat, index=rownames,
                       columns=["Estimate", "SE", "Z-stat", "P-val"])
  coefs["2.5_ci"]   = coefs["Estimate"] - 1.96 * coefs["SE"]
  coefs["97.5_ci"]  = coefs["Estimate"] + 1.96 * coefs["SE"]
  coefs["OR"]            = np.exp(coefs["Estimate"])
  coefs["OR_2.5_ci"]     = np.exp(coefs["2.5_ci"])
  coefs["OR_97.5_ci"]    = np.exp(coefs["97.5_ci"])
  coefs["Prob"]          = coefs["OR"]         / (1 + coefs["OR"])
  coefs["Prob_2.5_ci"]   = coefs["OR_2.5_ci"]  / (1 + coefs["OR_2.5_ci"])
  coefs["Prob_97.5_ci"]  = coefs["OR_97.5_ci"] / (1 + coefs["OR_97.5_ci"])
  coefs["Sig"] = pd.cut(coefs["P-val"],
                        bins=[-np.inf, 0.001, 0.01, 0.05, 0.1, np.inf],
                        labels=["***", "**", "*", ".", ""])
  ranef_var = vc.set_index("grp")[["vcov", "sdcor"]].rename(
              columns={"vcov": "Var", "sdcor": "Std"})
  # Match pymer4 convention: unnamed index so model_correct.py's
  # reset_index().rename({"index": "term"}) produces the term column.
  ranef_var.index.name = None
  ranef_var["Name"] = "(Intercept)"
  result = SimpleNamespace()
  result.coefs = coefs
  result.ranef_var = ranef_var
  return result

def plot_coefficients(coef_df, formula="predicted_doublet", feature_type="binary", cell_type="all"):
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
    # check if error is 
    for bar, (_,row) in zip(ax.patches, coef_df.iterrows()):
        bar_center = bar.get_y() + bar.get_height() / 2
        ci_low = row["2.5_ci"]
        ci_high = row["97.5_ci"]
        estimate = row["Estimate"]
    # check if error is greater than 1000
        if ci_low < -100 or ci_high > 100:
            print(f"Skipping error bar for {row['term']} due to extreme values: {ci_low}, {ci_high}")
            continue
 
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
    #formula_str = formula.replace("_", " ")
    plt.title(f"Correct ~ {feature_type} coefficients for {cell_type}")
    plt.xlabel("Log Odds")
    plt.ylabel("Feature")
    plt.tight_layout()
    plt.savefig(os.path.join(formula, f"{cell_type}_{feature_type}_coefficients.png"))
    
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



def build_formulas(predicted_meta, binary_features=None, 
                   continuous_features=None, 
                   doublet_features=["predicted_doublet", "doublet_score"]):
  
  formulas = {}
  for doublet_col in doublet_features:      # set type to binary or continuous
    formulas[doublet_col] = {
      "binary": None,
      "continuous": None
    } 
    for feature_type, feature_set in [("binary", binary_features), ("continuous", continuous_features)]:
      feature_cols = feature_set + [doublet_col]
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
        formula = " + ".join(valid_feature_cols)
        formulas[doublet_col][feature_type] = formula
      else:
        print(f"No valid features found for formula with {doublet_col}.")
  
  if not formulas:
    Warning("No valid formulas could be constructed.")
  return formulas