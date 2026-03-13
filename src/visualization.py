"""
Visualization utilities for flotation survival analysis.

This module provides functions for plotting model results:
  - Feature impact forest plots
  - Baseline cumulative hazard curves
  - Faceted boxplots of predicted floatability by stream and mineral
  - R² goodness-of-fit histograms
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_feature_impact(cph_model, figsize=(10, 8)):
    """
    Plot the impact of features on the flotation rate (hazard).

    Uses lifelines' built-in `plot()` method which generates a forest-style
    plot of coefficients with confidence intervals.

    Args:
        cph_model (lifelines.CoxPHFitter): Fitted Cox PH model.
        figsize (tuple): Figure size. Defaults to (10, 8).
    """
    plt.figure(figsize=figsize)
    cph_model.plot()
    plt.title('Feature Impact on Flotation Rate (Hazard)', fontsize=14)
    plt.xlabel('Coefficient (log-hazard ratio)')
    plt.tight_layout()
    plt.show()


def plot_baseline_hazard(cph_model, figsize=(10, 6)):
    """
    Plot the baseline cumulative hazard function.

    The baseline hazard represents the average flotation kinetics when
    all covariates are at their reference values (zero after centering).

    Args:
        cph_model (lifelines.CoxPHFitter): Fitted Cox PH model.
        figsize (tuple): Figure size. Defaults to (10, 6).
    """
    plt.figure(figsize=figsize)
    cph_model.baseline_cumulative_hazard_.plot()
    plt.title('Baseline Cumulative Hazard (Kinetic Rate Proxy)', fontsize=14)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Cumulative Hazard')
    plt.tight_layout()
    plt.show()


def get_stream(cls_name):
    """
    Extract the stream identifier from a class label.

    Maps class labels (e.g., CA20, CB50) to their stream code
    (CA, CB, CC, CD, TD) based on substring matching.

    Args:
        cls_name: Class label string.

    Returns:
        str: Stream code ('CA', 'CB', 'CC', 'CD', 'TD', or 'Other').
    """
    cls = str(cls_name).upper()
    if 'CA' in cls:
        return 'CA'
    if 'CB' in cls:
        return 'CB'
    if 'CC' in cls:
        return 'CC'
    if 'CD' in cls:
        return 'CD'
    if 'TD' in cls or 'TAIL' in cls:
        return 'TD'
    return 'Other'


def plot_floatability_boxplot(df, cph_model, data_model, figsize=None):
    """
    Create a faceted boxplot of predicted log-hazard by stream and
    dominant mineral.

    For each particle, the dominant mineral is determined by the mineral
    with the highest mass fraction. The predicted log partial hazard from
    the Cox model serves as a floatability risk score.

    Args:
        df (pd.DataFrame): Full DataFrame with 'Class' and mineral columns.
        cph_model (lifelines.CoxPHFitter): Fitted Cox PH model.
        data_model (pd.DataFrame): Model input data used for prediction.
        figsize (tuple, optional): Figure size override.
    """
    df = df.copy()
    df['Predicted_LogHazard'] = cph_model.predict_log_partial_hazard(data_model)

    df['Stream'] = df['Class'].apply(get_stream)

    mineral_columns = [
        'Apatite', 'Calcite', 'Dolomite', 'Albite', 'Quartz',
        'Actinolite', 'Chalcopyrite', 'Pyrite', 'Biotite',
    ]
    valid_minerals = [m for m in mineral_columns if m in df.columns]

    if valid_minerals:
        df['Dominant_Mineral'] = df[valid_minerals].idxmax(axis=1)
    else:
        df['Dominant_Mineral'] = 'Unknown'

    g = sns.catplot(
        data=df,
        x='Stream',
        y='Predicted_LogHazard',
        col='Dominant_Mineral',
        kind='box',
        col_wrap=3,
        height=4,
        aspect=1.2,
        order=['CA', 'CB', 'CC', 'CD', 'TD'],
        palette='Greys'
    )

    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(
        'Predicted Floatability (Log Hazard) by Stream and Mineral',
        fontsize=16
    )
    g.set_axis_labels("Stream", "Predicted Log Hazard (ln)")

    for ax in g.axes.flat:
        ax.axhline(0, ls='--', c='gray', alpha=0.5)

    plt.show()


def plot_r2_histogram(r2_values, figsize=(8, 5)):
    """
    Plot a histogram of R² values from individual particle kinetic fits.

    The R² distribution shows how well the survival model reproduces
    individual particle-level recovery kinetics.

    Args:
        r2_values (np.ndarray): Array of R² values, one per particle.
        figsize (tuple): Figure size. Defaults to (8, 5).
    """
    plt.figure(figsize=figsize)
    plt.hist(r2_values, bins=100, color='#555555', edgecolor='none')
    plt.title('$R^2$ Histogram (Vectorized)', fontsize=14)
    plt.xlim(0.8, 1.01)
    plt.xlabel('$R^2$')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.show()

    print(f"Mean R2: {np.nanmean(r2_values):.5f}")
