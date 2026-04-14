# normalize_impute.py
# Pipeline: log2 transform → median normalization → MinProb imputation
# Input:  PSME1_vs_Scram.csv
# Output: PSME1_vs_Scram_normalized.csv

import pandas as pd
import numpy as np
from scipy import stats

# 1. Load
df = pd.read_csv("PSME1_vs_Scram.csv", index_col="ProteinID")
sample_cols = df.columns.tolist()

print(f"[1/5] Loaded: {df.shape[0]} proteins x {df.shape[1]} samples")
print(f"      Missing values before processing: {df.isna().sum().sum()}")

# 2. Log2 transform
# Zero intensities represent non-detections, not true zeros — replace before log transform
df = df.replace(0, np.nan)
df_log = np.log2(df)

print("\n[2/5] Log2 transform complete")
print(df_log.describe().round(2))

# 3. Median normalization
# Shift each sample to a common median to correct for loading/technical variation
print("\n[3/5] Applying median normalization...")
sample_medians = df_log.median(axis=0)
global_median  = sample_medians.median()
df_norm        = df_log - sample_medians + global_median

print("      Sample medians post-normalization (should be equal):")
print(df_norm.median(axis=0).round(4))

# 4. MinProb imputation
# Missing values in LFQ data are typically MNAR (missing not at random) —
# proteins below detection threshold. MinProb models this by imputing from
# a downshifted normal distribution at the low end of each sample's intensity range.
# width_factor scales the imputation SD; shift_factor controls how far down to center it.
print("\n[4/5] Running MinProb imputation (MNAR model)...")

np.random.seed(42)

def minprob_impute(series, width_factor=0.3, shift_factor=1.8):
    observed  = series.dropna()
    mu        = observed.mean()
    sigma     = observed.std()
    # Center the imputation distribution below the observed signal
    center    = mu - shift_factor * sigma
    scale     = width_factor * sigma
    n_missing = series.isna().sum()
    imputed   = series.copy()
    # Draw replacement values from the downshifted Gaussian
    imputed[series.isna()] = np.random.normal(center, scale, n_missing)
    return imputed

df_imputed = df_norm.apply(minprob_impute, axis=0)

print(f"      Missing values after imputation: {df_imputed.isna().sum().sum()}")
print("\n      Final matrix summary:")
print(df_imputed.describe().round(2))

# 5. Save outputs
print("\n[5/5] Saving outputs...")

df_imputed.to_csv("PSME1_vs_Scram_normalized.csv")
print("      Saved: PSME1_vs_Scram_normalized.csv")

# Per-sample QC stats for downstream review
qc = pd.DataFrame({
    "raw_median_log2":      df_log.median(),
    "norm_median_log2":     df_norm.median(),
    "n_missing_pre_impute": df_norm.isna().sum(),
    "mean_final":           df_imputed.mean(),
    "sd_final":             df_imputed.std(),
})
qc.to_csv("QC_sample_stats.csv")
print("      Saved: QC_sample_stats.csv")
print("\nNormalization pipeline complete.")