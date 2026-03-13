# Froth Flotation Survival Analysis

A Cox Proportional Hazards model with bootstrap LASSO regularization applied to mineral particle flotation data. This project models individual particle recovery kinetics by framing flotation as a survival problem: concentrates are "events" and tailings are "censored" observations.

## Methodology

| Step | Technique | Purpose |
|------|-----------|---------|
| 1. Data Balancing | Class-stratified downsampling | Equal representation (~60k per stream) |
| 2. Survival Mapping | Time/event encoding | CA→0.75min, CB→1.5min, CC→3min, CD→6min, TD→6min (censored) |
| 3. Feature Selection | VIF-based elimination | Remove multicollinear features (threshold VIF > 10) |
| 4. Model Fitting | Bootstrap Cox PH + LASSO | 1000 iterations × 50k subsamples, penalizer=0.02, L1 ratio=1.0 |
| 5. Coefficient Estimation | Bootstrap aggregation | Mean coefficients with 95% CIs from bootstrap distribution |

### Key Concepts

- **Cox Proportional Hazards**: Models the hazard (instantaneous flotation rate) as a function of particle properties. A positive coefficient means faster flotation; negative means slower.
- **LASSO Regularization (L1)**: Drives unimportant feature coefficients to exactly zero, performing automatic feature selection.
- **Bootstrap Aggregation**: Fits the model 1000 times on resampled data to estimate coefficient stability and confidence intervals.
- **VIF (Variance Inflation Factor)**: Pre-filters features to remove redundant/collinear variables before modeling.

## Project Structure

```
├── Data_Analysis_Project.ipynb   # Main Jupyter notebook (original analysis)
├── src/
│   ├── __init__.py               # Package init
│   ├── data_loader.py            # Data loading, class balancing, survival mapping
│   ├── feature_engineering.py    # Feature selection, VIF-based collinearity removal
│   ├── model.py                  # Bootstrap Cox PH fitting, result saving
│   └── visualization.py         # Plotting functions (forest, hazard, boxplots, R²)
├── RESULTS.md                    # Detailed explanation of findings
├── requirements.txt              # Python dependencies
├── LICENSE                       # MIT License
├── .gitignore                    # Excludes data files, caches, etc.
└── README.md                     # This file
```

## Data

The project uses particle-level mineral composition and shape data from flotation experiments. Data files are **not included** in this repository due to their size (~700MB total):

- `Traindata_20um.csv` (~207MB) — Primary dataset used in the analysis
- `Traindata.csv` (~485MB) — Full dataset

Place data files in the project root directory before running the notebook.

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd "Data Analysis 2"

# Install dependencies
pip install -r requirements.txt
```

### Requirements

- Python 3.10+
- lifelines ≥ 0.30.0
- pandas ≥ 2.1
- numpy ≥ 1.14.0
- matplotlib ≥ 3.0
- seaborn ≥ 0.12.0
- joblib ≥ 1.0.0
- statsmodels ≥ 0.13.0
- scikit-learn ≥ 1.0.0

## Usage

### Jupyter Notebook

Open `Data_Analysis_Project.ipynb` in Jupyter or Google Colab and run cells sequentially.

### Python API

```python
from src.data_loader import load_data, balance_classes, map_survival_variables
from src.feature_engineering import select_features, remove_collinear
from src.model import fit_bootstrap_cox, fit_final_model, save_bootstrap_results
from src.visualization import (
    plot_feature_impact,
    plot_baseline_hazard,
    plot_floatability_boxplot,
    plot_r2_histogram
)

# 1. Load and prepare data
df = load_data('Traindata_20um.csv')
df = balance_classes(df, target_total_n=300000)
df = map_survival_variables(df)

# 2. Feature engineering
feature_cols, data_model = select_features(df)
model_data = remove_collinear(data_model)

# 3. Bootstrap model fitting
boot_results = fit_bootstrap_cox(model_data, n_bootstrap=1000)
save_bootstrap_results(boot_results)
cph_final = fit_final_model(model_data, boot_results)

# 4. Visualize results
plot_feature_impact(cph_final)
plot_baseline_hazard(cph_final)
plot_floatability_boxplot(df, cph_final, model_data)
```

## Results Summary

See [RESULTS.md](RESULTS.md) for detailed findings. Key highlights:

| Feature | Coefficient | Hazard Ratio | Interpretation |
|---------|------------|--------------|----------------|
| **Apatite** | +0.55 | 1.74 | 74% faster flotation (ore mineral — most floatable) |
| **Phlogopite.surf** | −0.99 | 0.37 | 63% slower flotation (hydrophilic surface suppresses recovery) |
| **Quartz.surf** | −0.52 | 0.59 | 41% slower flotation |
| **Albite.surf** | −0.52 | 0.60 | 40% slower flotation |
| **Orthoclase.surf** | −0.30 | 0.74 | 26% slower flotation |
| **ECD** | +0.02 | 1.02 | Larger particles float slightly faster |

- **Model concordance**: 0.62
- **Mean R²** (individual kinetics): 0.96
- **Successful bootstrap fits**: 1000/1000

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
