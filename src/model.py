"""
Cox Proportional Hazards model fitting with bootstrap LASSO.

This module implements parallelized bootstrap fitting of a Cox PH model
with L1 (LASSO) regularization. The bootstrap estimates are aggregated
to produce stable coefficient estimates and confidence intervals.
"""

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from lifelines import CoxPHFitter


def fit_one_bootstrap(model_data, all_features, seed,
                      subsample_n=50000, penalizer=0.02, l1_ratio=1.0):
    """
    Fit a single bootstrap iteration of the Cox PH model.

    Draws a bootstrap sample (with replacement), checks for zero-variance
    features in that sample, and fits a LASSO-regularized Cox PH model.

    Args:
        model_data (pd.DataFrame): The full modeling dataset with features + 'T' + 'E'.
        all_features (list): List of all feature column names.
        seed (int): Random seed for this bootstrap sample.
        subsample_n (int): Maximum rows per bootstrap sample. Defaults to 50000.
        penalizer (float): Regularization strength. Defaults to 0.02.
        l1_ratio (float): L1/L2 mixing ratio (1.0 = pure LASSO). Defaults to 1.0.

    Returns:
        pd.Series or None: Fitted coefficients reindexed to all_features
                           (missing features filled with 0), or None on failure.
    """
    boot_sample = model_data.sample(
        n=min(subsample_n, len(model_data)),
        replace=True,
        random_state=seed
    )

    # Per-sample variance check
    sample_stds = boot_sample[all_features].std()
    valid_features = sample_stds[sample_stds > 1e-6].index.tolist()
    boot_clean = boot_sample[valid_features + ['T', 'E']]

    try:
        cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
        cph.fit(boot_clean, duration_col='T', event_col='E')
        return cph.params_.reindex(all_features, fill_value=0.0)
    except Exception:
        return None


def fit_bootstrap_cox(model_data, n_bootstrap=1000, subsample_n=50000,
                      penalizer=0.02, l1_ratio=1.0, n_jobs=-1, verbose=5):
    """
    Run parallelized bootstrap fitting of a Cox PH model with LASSO.

    Each iteration draws a bootstrap subsample, fits the model, and
    collects coefficient estimates. Results are aggregated to compute
    mean coefficients, standard errors, and 95% confidence intervals.

    Args:
        model_data (pd.DataFrame): Modeling dataset with features + 'T' + 'E'.
        n_bootstrap (int): Number of bootstrap iterations. Defaults to 1000.
        subsample_n (int): Rows per bootstrap sample. Defaults to 50000.
        penalizer (float): Regularization strength. Defaults to 0.02.
        l1_ratio (float): L1/L2 mixing ratio. Defaults to 1.0.
        n_jobs (int): Number of parallel jobs (-1 = all CPUs). Defaults to -1.
        verbose (int): Verbosity level for joblib. Defaults to 5.

    Returns:
        pd.DataFrame: DataFrame of bootstrap coefficient estimates,
                      one row per successful iteration.
    """
    all_features = [c for c in model_data.columns if c not in ['T', 'E']]

    print(f"Starting parallel bootstrap ({n_bootstrap} iters, "
          f"{subsample_n} rows each)...")

    results = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(fit_one_bootstrap)(
            model_data, all_features, i, subsample_n, penalizer, l1_ratio
        )
        for i in range(n_bootstrap)
    )

    bootstrapped_betas = [r for r in results if r is not None]
    boot_results = pd.DataFrame(bootstrapped_betas)

    print(f"\nSuccessful fits: {len(boot_results)}/{n_bootstrap}")
    print("\n=== Bootstrap + LASSO Results ===")
    print("\nMean Coefficients:")
    print(boot_results.mean())
    print("\n95% Confidence Intervals:")
    print(boot_results.quantile([0.025, 0.975]).T)
    print("\nStandard Errors:")
    print(boot_results.std())

    return boot_results


def save_bootstrap_results(boot_results, prefix='bootstrap'):
    """
    Save bootstrap results to CSV files.

    Generates four files:
      - {prefix}_mean_coefficients.csv
      - {prefix}_std_errors.csv
      - {prefix}_confidence_intervals.csv
      - {prefix}_summary.csv

    Args:
        boot_results (pd.DataFrame): Bootstrap coefficient estimates.
        prefix (str): Filename prefix. Defaults to 'bootstrap'.
    """
    boot_results.mean().to_csv(
        f'{prefix}_mean_coefficients.csv', header=['mean_coef']
    )
    boot_results.std().to_csv(
        f'{prefix}_std_errors.csv', header=['std_error']
    )
    boot_results.quantile([0.025, 0.975]).T.to_csv(
        f'{prefix}_confidence_intervals.csv'
    )

    summary = pd.DataFrame({
        'Mean Coef': boot_results.mean(),
        'Std Error': boot_results.std(),
        'CI 2.5%': boot_results.quantile(0.025),
        'CI 97.5%': boot_results.quantile(0.975)
    })
    summary.to_csv(f'{prefix}_summary.csv')
    print(f"Saved to {prefix}_summary.csv")
    print(summary)


def fit_final_model(model_data, boot_results, penalizer=0.02, l1_ratio=1.0):
    """
    Fit the final Cox PH model on the full dataset and override
    its coefficients with the bootstrap mean estimates.

    This approach combines the stability of bootstrap estimates with
    the baseline hazard estimated from the full dataset.

    Args:
        model_data (pd.DataFrame): Full modeling dataset.
        boot_results (pd.DataFrame): Bootstrap coefficient estimates.
        penalizer (float): Regularization strength. Defaults to 0.02.
        l1_ratio (float): L1/L2 mixing ratio. Defaults to 1.0.

    Returns:
        lifelines.CoxPHFitter: Fitted model with bootstrap-averaged coefficients.
    """
    all_features = [c for c in model_data.columns if c not in ['T', 'E']]

    # Variance check on full data
    f_stds = model_data[all_features].std()
    valid_final = f_stds[f_stds > 1e-6].index.tolist()
    final_data = model_data[valid_final + ['T', 'E']]

    # Fit final model on all rows
    cph_final = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
    cph_final.fit(final_data, duration_col='T', event_col='E')

    # Override coefficients with bootstrap means
    cph_final.params_ = boot_results.mean().reindex(
        cph_final.params_.index, fill_value=0.0
    )

    print("\n=== Final Model Summary (using Bootstrap means) ===")
    cph_final.print_summary()

    return cph_final
