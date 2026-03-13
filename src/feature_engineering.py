"""
Feature engineering and selection for flotation survival analysis.

This module provides utilities for selecting numeric features,
removing constant columns, and eliminating multicollinear features
using Variance Inflation Factor (VIF) analysis.
"""

import pandas as pd
import numpy as np
from statsmodels.stats.outliers_influence import variance_inflation_factor


# Columns that should never be used as features
EXCLUDE_COLS = ['Class', 'T', 'E', 'sample_weights']


def select_features(df, exclude_cols=None):
    """
    Identify numeric feature columns and drop any that are constant
    (zero variance).

    Args:
        df (pd.DataFrame): Input DataFrame.
        exclude_cols (list, optional): Column names to exclude.
            Defaults to EXCLUDE_COLS.

    Returns:
        tuple: (feature_cols, data_model) where feature_cols is a list
               of selected feature names and data_model is the filtered
               DataFrame with features + ['T', 'E'], NaN rows dropped.
    """
    if exclude_cols is None:
        exclude_cols = EXCLUDE_COLS

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    feature_cols = [c for c in numeric_cols if c not in exclude_cols]

    # Drop truly constant columns (std == 0)
    stds = df[feature_cols].std()
    constant_cols = stds[stds == 0].index.tolist()
    if constant_cols:
        print(f"Dropping {len(constant_cols)} constant columns: {constant_cols}")
        feature_cols = [c for c in feature_cols if c not in constant_cols]

    data_model = df[feature_cols + ['T', 'E']].dropna()
    print(f"\nUsing {len(feature_cols)} features on {len(data_model)} rows")
    print(feature_cols)
    return feature_cols, data_model


def remove_collinear(df_in, threshold=10.0, sample_n=50000):
    """
    Iteratively remove features with high Variance Inflation Factor (VIF)
    to address multicollinearity.

    At each step, the feature with the highest VIF above the threshold
    is dropped, until all remaining features have VIF ≤ threshold.

    A subsample of the data is used for VIF computation to improve speed.

    Args:
        df_in (pd.DataFrame): DataFrame with feature columns + 'T' and 'E'.
        threshold (float): VIF threshold. Features with VIF > threshold
            are removed. Defaults to 10.0.
        sample_n (int): Number of rows to subsample for VIF calculation.
            Defaults to 50000.

    Returns:
        pd.DataFrame: Filtered DataFrame with collinear features removed,
                       retaining 'T' and 'E' columns.
    """
    df_s = df_in.sample(
        n=min(sample_n, len(df_in)),
        random_state=42
    )
    X = df_s.drop(['T', 'E'], axis=1)

    while True:
        vifs = []
        for i in range(len(X.columns)):
            try:
                vifs.append(variance_inflation_factor(X.values, i))
            except Exception:
                vifs.append(float('inf'))

        v = pd.DataFrame({'f': X.columns, 'V': vifs})
        mx = v['V'].max()

        if pd.isna(mx) or mx == float('inf') or mx > threshold:
            feat = v.sort_values(
                'V', ascending=False, na_position='first'
            )['f'].iloc[0]
            X = X.drop(columns=[feat])
            print(f"  Dropped collinear: {feat}")
        else:
            break

    return df_in[X.columns.tolist() + ['T', 'E']]
