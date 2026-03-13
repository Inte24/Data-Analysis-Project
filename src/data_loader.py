"""
Data loading and preprocessing utilities for flotation survival analysis.

This module handles loading CSV data, balancing class distributions,
and mapping flotation class labels to survival analysis variables
(time T and event E).
"""

import pandas as pd
import numpy as np


def load_data(filename):
    """
    Load particle data from a CSV file.

    Args:
        filename (str): Path to the CSV file containing particle data.

    Returns:
        pd.DataFrame: The loaded DataFrame.
    """
    print(f"Loading {filename}...")
    df = pd.read_csv(filename)
    print(f"Original shape: {df.shape}")
    return df


def balance_classes(df, target_total_n=300000, random_state=42):
    """
    Balance the dataset by downsampling larger classes to achieve
    approximately equal representation.

    Each class is sampled to `target_total_n // n_classes` rows.
    Classes with fewer rows than the target are kept as-is.

    Args:
        df (pd.DataFrame): Input DataFrame with a 'Class' column.
        target_total_n (int): Desired total number of samples. Defaults to 300000.
        random_state (int): Random seed for reproducibility. Defaults to 42.

    Returns:
        pd.DataFrame: Balanced DataFrame with reset index.
    """
    relevant_classes = [c for c in df['Class'].unique() if isinstance(c, str)]
    n_classes = len(relevant_classes)
    samples_per_class = target_total_n // n_classes

    print(f"\nBalancing: Aiming for ~{samples_per_class} samples per class "
          f"for {n_classes} classes.")

    balanced_dfs = []
    for cls in relevant_classes:
        cls_df = df[df['Class'] == cls]
        if len(cls_df) >= samples_per_class:
            balanced_dfs.append(
                cls_df.sample(n=samples_per_class, random_state=random_state)
            )
        else:
            balanced_dfs.append(cls_df)

    result = pd.concat(balanced_dfs).reset_index(drop=True)
    print(f"Balanced shape: {result.shape}")
    print("\nClass distribution in balanced set:")
    print(result['Class'].value_counts())
    return result


def get_time_event(cls_name):
    """
    Map a flotation class label to survival analysis variables.

    Concentrate streams (CA, CB, CC, CD) are treated as events (E=1)
    with cumulative flotation times. Tailings (TD/TAIL) are right-censored
    observations (E=0) at the final collection time.

    Time mapping:
        - CA → T=0.75 min (E=1)
        - CB → T=1.50 min (E=1)
        - CC → T=3.00 min (E=1)
        - CD → T=6.00 min (E=1)
        - TD/TAIL → T=6.00 min (E=0, censored)

    Args:
        cls_name: The class label (e.g., 'CA20', 'CB50', 'TD20').

    Returns:
        tuple: (time, event) pair, or (None, None) if class is unknown.
    """
    s = str(cls_name).upper()

    if 'CA' in s:
        return 0.75, 1
    elif 'CB' in s:
        return 1.50, 1
    elif 'CC' in s:
        return 3.00, 1
    elif 'CD' in s:
        return 6.00, 1
    elif 'TD' in s or 'TAIL' in s:
        return 6.00, 0  # Censored at final collection time
    else:
        return None, None


def map_survival_variables(df):
    """
    Add survival analysis columns (T, E) to the DataFrame based on
    the 'Class' column.

    T (Time): Cumulative flotation time in minutes.
    E (Event): 1 = floated (recovered in concentrate), 0 = not floated (censored).

    Args:
        df (pd.DataFrame): DataFrame with a 'Class' column.

    Returns:
        pd.DataFrame: DataFrame with added 'T' and 'E' columns,
                       rows with unknown classes removed.
    """
    df[['T', 'E']] = df['Class'].apply(
        lambda x: pd.Series(get_time_event(x))
    )
    df = df.dropna(subset=['T', 'E'])
    return df
