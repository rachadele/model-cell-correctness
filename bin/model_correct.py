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
from collections import defaultdict
from pymer4.models import Lmer


def parse_arguments():
    parser = argparse.ArgumentParser(description="Correct model predictions based on reference keys and mapping file.")
    parser.add_argument('--predicted_meta', type=str, help="Path to the predicted metadata file", default="/space/grp/rschwartz/rschwartz/model_cell_correctness/work/5a/0c743fb7e8606b0439689f4a940e57/Immune_predicted_meta_subset.tsv")
    parser.add_argument('--mapping_file', type=str, help="Path to the cell type hierarchy file", default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_keys', type=str, nargs="+", help="Reference keys to map", default=["subclass", "class", "family"])
    parser.add_argument('--cell_type', type=str, help="Cell type to analyze", default="Immune")
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
    binary_features = ["ribo_outlier", "mito_outlier", "hb_outlier", "counts_outlier", "genes_outlier", "umi_outlier"]
    continuous_features = ["pct_counts_ribo", "pct_counts_mito", "pct_counts_hb", "log1p_n_genes_by_counts", "log1p_total_counts"]

    for col in binary_features + continuous_features + ["predicted_doublet", "doublet_score"]:
        if pd.api.types.is_bool_dtype(predicted_meta[col]):
            predicted_meta[col] = predicted_meta[col].astype(int)
        elif pd.api.types.is_numeric_dtype(predicted_meta[col]):
            predicted_meta[col] = predicted_meta[col].astype(float)
            
            
    # Check whether to include random effect
    random_effect_column = "study"
    use_random_effect = predicted_meta[random_effect_column].nunique() > 1

    
    outcomes = [f"correct_{key}" for key in ref_keys]
    columns_to_keep = binary_features + continuous_features + outcomes + ["study"] + ["sample_id"] + ["predicted_doublet", "doublet_score"]
    predicted_meta = predicted_meta[columns_to_keep]
    
    
    # Build formulas (without adding random effect here)
    formulas = build_formulas(predicted_meta, 
                              binary_features=binary_features, 
                              continuous_features=continuous_features,
                              doublet_features=["predicted_doublet"])


    models = {}
    for formula_name, formula_dict in formulas.items():
        models[formula_name] = {}
        for feature_type, formula in formula_dict.items():
            if formula is None:
                continue
            models[formula_name][feature_type] = {}
            
            for key in ref_keys:
                true_col = f"{key}"
                outcome = f"correct_{true_col}"
                predicted_meta[outcome] = predicted_meta[outcome].astype(int)
                new_formula = f"{outcome} ~ {formula}"
                models[formula_name][feature_type][key] = {}
                try:
                    if use_random_effect:
                        full_formula = f"{new_formula} + (1|{random_effect_column})"
                        model = run_llmer(predicted_meta, formula=full_formula, family="binomial")
                    else:
                        full_formula = new_formula
                        model = smf.logit(formula=full_formula, data=predicted_meta).fit(disp=False)
                    
                    models[formula_name][feature_type][key]["fit"] = model
                    models[formula_name][feature_type][key]["formula"] = full_formula

                except Exception as e:
                    print(f"Error fitting model for {key} with formula {full_formula}: {e}")
                    continue

    for formula_name, feature_type in models.items():
        for feature_type, key_model_dict in feature_type.items():
            coef_list = []

            for key, subdict in key_model_dict.items():
                if "fit" not in subdict:
                    print(f"Skipping {key} in {formula_name} as fit is not available.")
                    continue
                # doesn't work when entire coef list is empty

                model = subdict["fit"]

                if hasattr(model, "coefs"):  # pymer4
                    model_df = model.coefs.reset_index().rename(columns={"index": "term"})
                else:  # statsmodels
                    sm_df = model.summary2().tables[1].copy()
                    sm_df = sm_df.reset_index().rename(columns={
                        "index": "term", 
                        "Coef.": "Estimate", 
                        "[0.025": "2.5_ci", 
                        "0.975]": "97.5_ci"
                    })
                    model_df = sm_df[["term", "Estimate", "2.5_ci", "97.5_ci"]]

                for _, row in model_df.iterrows():
                    row_dict = row.to_dict()
                    row_dict["key"] = key
                    row_dict["formula"] = subdict["formula"]
                    coef_list.append(row_dict)

                # Handle random effects variance from pymer4
                if use_random_effect and hasattr(model, "ranef_var"):
                    ranef_var = model.ranef_var.reset_index().rename(columns={"index": "term"})
                    for _, row in ranef_var.iterrows():
                        row_dict = row.to_dict()
                        row_dict["Estimate"] = row_dict["Var"]
                        row_dict["2.5_ci"] = row_dict["Var"] - 2 * row_dict["Std"]
                        row_dict["97.5_ci"] = row_dict["Var"] + 2 * row_dict["Std"]
                        row_dict["key"] = key
                        row_dict["formula"] = subdict["formula"]
                        coef_list.append(row_dict)

            # Save and plot
            if not coef_list:
                print(f"No coefficients found for {formula_name} with feature type {feature_type}.")
                continue
            coef_df = pd.DataFrame(coef_list)
            os.makedirs(formula_name, exist_ok=True)
            coef_df.to_csv(os.path.join(formula_name,f"{cell_type}_{feature_type}_coefficients.tsv"), sep="\t", index=False)
            plot_coefficients(coef_df, formula=formula_name, cell_type=cell_type, feature_type=feature_type)

if __name__ == "__main__":
    main()
