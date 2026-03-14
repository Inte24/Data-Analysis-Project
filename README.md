# Froth Flotation Survival Analysis

A Cox Proportional Hazards model with bootstrap LASSO regularization applied to mineral particle flotation data. This project models individual particle recovery kinetics by framing flotation as a survival problem: concentrates are "events" and tailings are "censored" observations.

## Understanding the Data

Each row in the dataset represents a single mineral particle analyzed by automated mineralogy. For each particle, the dataset contains:

- **Mineral mass fractions**: Bulk composition (e.g., Apatite, Quartz, Calcite, Pyrite, etc.)
- **Mineral surface fractions** (`.surf`): Proportion of the particle's exposed surface covered by each mineral, this is what determines flotation behavior
- **Shape descriptors**: AspectRatio, Solidity, ECD (Equivalent Circle Diameter in µm)
- **Class label**: Encodes both the flotation stream and size fraction

### Class Labels

The class label combines a stream code with a size fraction. For example, `CA20` = a particle of the −20 µm size fraction that floated in 0.75 min.

| Stream | Cumulative Time | Meaning |
|--------|----------------|---------|
| CA | 0.75 min | First concentrate: fastest floating particles |
| CB | 1.50 min | Second concentrate |
| CC | 3.00 min | Third concentrate |
| CD | 6.00 min | Fourth concentrate: slowest floating particles |
| TD | Never floated | Tailings: reported to waste (right-censored) |

| Size Fraction | Particle Size |
|---------------|---------------|
| 20 | −20 µm |
| 2032 | +20 to −32 µm |
| 3250 | +32 to −50 µm |
| 50 | +50 µm |

> This analysis uses only the −20 µm fraction (`Traindata_20um.csv`) for computational efficiency.

### Data Files

Data files are not included in this repository due to their size (~700MB total):

- `Traindata_20um.csv` (~207MB) — Primary dataset (796,130 particles, 54 columns) - https://drive.google.com/file/d/1GVB0wzckLw6iFmtLfhzQ24UNjgx6fzV4/view?usp=drive_link
- `Traindata.csv` (~485MB) — Full dataset with all size fractions - https://drive.google.com/file/d/1JaUAQ9oif9Ai7n_gdvylriAvjs8Ine_D/view?usp=drive_link

Place data files in the project root directory before running the notebook.

---

## Data Preprocessing

1. **Size fraction selection**: Only −20 µm particles are used, reducing complexity while maintaining statistical power (796,130 particles)
2. **Class balancing**: Downsample to 300,000 total rows with ~60,000 particles per class (CA20, CB20, CC20, CD20, TD20). This prevents the model from being biased toward overrepresented classes
3. **Survival variable mapping**: Map each class label to:
   - **T (Time)**: Cumulative flotation time when the particle was recovered
   - **E (Event)**: 1 = floated (recovered in concentrate/ Death), 0 = not floated (tailings, Survived)

---

## Mathematical Formulation

### From Flotation Kinetics to Survival Analysis

Classical first-order flotation kinetics describes recovery as:

$$F(t) = \frac{R(t)}{R_{max}} = 1 - e^{-kt}$$

where *k* is the kinetic rate constant. The survival function (probability of NOT being recovered by time *t*) is:

$$S(t) = 1 - F(t) = e^{-kt}$$

The hazard function (instantaneous flotation rate) is constant for first-order kinetics:

$$h(t) = -\frac{d}{dt} \ln S(t) = k$$

### Cox Proportional Hazards Extension

The Cox PH model extends this by allowing the hazard to depend on particle properties x = (mineralogy, surface, size):

$$h(t \mid \mathbf{x}) = h_0(t) \cdot \exp\left(\sum_{i=1}^{n} \beta_i x_i\right)$$

where:
- *h₀(t)* is the baseline hazard (average flotation kinetics)
- *βᵢ* are the coefficients — positive = faster flotation, negative = slower
- *xᵢ* are the particle properties (mineral fractions, surface fractions, ECD)

The corresponding survival function becomes:

$$S(t \mid \mathbf{x}) = \exp\left(-H_0(t) \cdot e^{\sum \beta_i x_i}\right)$$

where *H₀(t)* is the baseline cumulative hazard. This framework allows us to model individual particle kinetics based on their measured properties.

---

## Data Cleaning & Feature Selection

The Cox PH model is sensitive to multicollinearity: correlated features cause unstable coefficient estimates and singular matrix errors. We address this with a three-step pipeline:

### Step 1: Remove Constant Features
- Drop columns with zero standard deviation (e.g., Clinochlore, Clinochlore.surf, these minerals were absent in the −20 µm fraction)
- Constant values would cause a singular matrix in the Cox model's partial likelihood optimization

### Step 2: VIF-Based Collinearity Removal
- Variance Inflation Factor (VIF) measures how much a feature's variance is inflated by correlation with other features
- A VIF > 10 means that 90% of the feature's variance can be explained by other features — it is essentially a clone of another variable
- VIF is calculated on a 50,000-row subsample for computational efficiency (computing VIF on all 300k rows is expensive, though a correlation matrix is an alternative)
- Features are dropped iteratively in a loop, because removing one feature changes the VIF of all remaining features
- Result: 51 features → 24 features after VIF filtering

### Step 3: LASSO Feature Selection
- L1 regularization (penalizer=0.02, l1_ratio=1.0) drives unimportant coefficients to exactly zero
- Of the 24 features surviving VIF, LASSO further identifies ~6 significant predictors with non-negligible coefficients

---

## Combining Bootstrap and LASSO Regularization

Before bootstrapping, it was critical to find the correct penalizer (regularization strength) for the LASSO. The `l1_ratio` was fixed at 1.0 (pure LASSO), and six penalizer values were tested:

### Penalizer Comparison

| Penalizer | Non-zero Features | Behavior |
|-----------|------------------|----------|
| 0.001 | 9 | Too weak, retains noise features with small, unreliable coefficients |
| 0.005 | ~7 | Still retains some marginal features |
| 0.01 | ~6 | Getting closer, but some weak signals remain |
| **0.02** | **5** | **Sweet spot, retains only strong, interpretable predictors** |
| 0.05 | ~3 | Starting to over-penalize, loses some real signals |
| 0.1 | 2 | Too aggressive, nearly all coefficients driven to ≈0, only Quartz and ECD survive |

### Why Penalizer = 0.02?

**Too low (0.001)**: The model retains 9 features with non-zero coefficients, including noisy ones. For example at penalizer=0.001, Actinolite (−0.56), Solidity (−0.37), and Chalcopyrite (−0.008) all have non-zero coefficients, but these are not stable across bootstrap iterations, they are likely artifacts of overfitting.

**Too high (0.1)**: The model shrinks almost everything to zero. At penalizer=0.1, features like Albite (−3.4e-06) and Biotite (−3.5e-06) that should be meaningful predictors are effectively eliminated. Only Quartz.surf (−0.034) and ECD (0.012) survive, the model is too aggressive and loses real signal.

**Just right (0.02)**: The model retains 5 features with meaningful coefficients:

| Feature | Coefficient at λ=0.02 |
|---------|----------------------|
| Albite | −0.81 |
| Quartz | −0.81 |
| AspectRatio | +0.35 |
| ECD | +0.03 |
| Pyrite | −0.03 |

This produces a sparse, interpretable model where retained features have coefficients large enough to be physically meaningful, while noise features are cleanly driven to zero. The penalizer=0.02 is the sweet spot that balances sparsity (few features) with signal retention (keeps real predictors).

---

## Bootstrap Aggregation

### What is Bootstrapping?

Instead of fitting the Cox PH model once on the full dataset, we fit it 1000 times on random subsamples (50,000 rows each, drawn with replacement). Each iteration produces a set of coefficient estimates. The final coefficients are the mean across all 1000 fits. This also makes sure to mimic real life situations where not each class will be balanced and we might get an instant where one class outweighs the other. 

### Why Bootstrap?

1. **Coefficient stability**: A single model fit can be sensitive to the specific rows in the dataset. Bootstrapping reveals whether coefficients are consistently positive/negative across different data subsets, or if they fluctuate around zero (unreliable)
2. **Robust confidence intervals**: The 2.5th and 97.5th percentiles of the bootstrap distribution provide empirical 95% confidence intervals that do not rely on normality assumptions
3. **Variance estimation**: The standard deviation across bootstrap iterations gives a direct estimate of coefficient uncertainty
4. **LASSO stability selection**: With L1 regularization, some features may be selected in some bootstrap iterations and not others. Features that are consistently non-zero across many iterations are more reliable than those that appear sporadically

### Choosing Bootstrap Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Subsample size | 50,000 | Computational trade-off (see below) |
| Iterations | 1,000 | Statistical convergence (see below) |
| Sampling | With replacement | Standard bootstrap methodology |
| Parallelization | All CPU cores | `joblib.Parallel` with `n_jobs=-1` |

**Why n = 50,000 per iteration (not 300,000)?**

Fitting a Cox PH model requires optimizing a partial likelihood over all risk sets at each event time. The computational cost scales super-linearly with dataset size, fitting on 300,000 rows is not simply 6× slower than 50,000 rows, it is substantially more expensive due to the sorting and pairwise comparisons involved. Using 50,000-row subsamples:
- Keeps each individual fit fast (~seconds rather than minutes)
- Still provides enough statistical power to estimate coefficients reliably (50k rows with 24 features is well above the rule-of-thumb minimum of ~10 events per feature)
- Reduces the risk of convergence failures that can occur with LASSO on very large datasets
- Makes 1000 iterations feasible: 1000 fits on 50k rows is faster than 1000 fits on 300k rows, which would take days

**Why 1,000 iterations?**

- By the law of large numbers, the mean coefficient estimate converges to the true population value as the number of iterations increases
- The standard error of the bootstrap mean scales as σ/√n, with 1000 iterations, the uncertainty in the mean coefficient is ~3% of the standard deviation across iterations (1/√1000 ≈ 0.032)
- Beyond 1000, there are diminishing returns: going from 1000 to 10,000 iterations would only improve precision by a factor of √10 ≈ 3.2×, at 10× the computational cost
- 1000 iterations also gives reliable percentile-based confidence intervals, the 2.5th percentile is estimated from ~25 data points, which is sufficient for stable CI bounds

### Bootstrap Workflow

1. Draw 50,000 rows with replacement from the 300,000-row balanced dataset
2. Check for zero-variance features in this specific sample (rare features may have no variance in a subsample)
3. Fit Cox PH with LASSO on the subsample
4. Record the 24 coefficients
5. Repeat 1000 times (parallelized across all CPU cores)
6. Aggregate: take the mean of each coefficient across all iterations
7. Override the final model's coefficients with these bootstrap means

All 1000 iterations converged successfully. The narrow confidence intervals (particularly for Apatite, Phlogopite.surf, Quartz.surf) confirm these are robust, reliable predictors.

---

## Project Structure

```
├── Data_Analysis_Project.ipynb   # Main Jupyter notebook (original analysis)
├── Data_Analysis_Pipeline.ipynb  # Clean notebook using src/ pipeline
├── src/
│   ├── __init__.py               # Package init
│   ├── data_loader.py            # Data loading, class balancing, survival mapping
│   ├── feature_engineering.py    # Feature selection, VIF-based collinearity removal
│   ├── model.py                  # Bootstrap Cox PH fitting, result saving
│   └── visualization.py         # Plotting functions (forest, hazard, boxplots, R²)
├── figures/                      # Generated plots and visualizations
├── RESULTS.md                    # Detailed results with figure interpretations
├── requirements.txt              # Python dependencies
├── LICENSE                       # MIT License
├── .gitignore                    # Excludes data files, caches, etc.
└── README.md                     # This file
```

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
| **Apatite** | +0.55 | 1.74 | 74% faster flotation (ore mineral, most floatable) |
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
