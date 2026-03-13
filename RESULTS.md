# Results & Interpretation

## Model Overview

A Cox Proportional Hazards model with LASSO regularization was fitted to ~300,000 balanced mineral particle observations from froth flotation experiments. The model treats flotation as a survival process where:

- **Event** = particle is recovered into concentrate (streams CA, CB, CC, CD)
- **Censored** = particle reports to tailings (TD) and is not recovered

Bootstrap aggregation (1000 iterations, 50k subsamples each) was used to ensure stable coefficient estimates. All 1000 bootstrap fits converged successfully.

## Model Performance

| Metric | Value |
|--------|-------|
| Number of observations | 300,000 |
| Events observed | 240,000 |
| Right-censored | 60,000 |
| Concordance index | 0.62 |
| Partial AIC | 5,735,176.37 |
| Mean R² (particle kinetics) | 0.96 |

The concordance of 0.62 indicates the model correctly ranks particles by flotation speed 62% of the time — a reasonable result for individual particle data where stochastic effects are significant.

The mean R² of 0.96 demonstrates that the model reproduces individual particle recovery kinetics very well, validating the survival analysis framework for flotation.

## Significant Features

After LASSO feature selection, six features had non-negligible coefficients. Several features (Chalcopyrite, Diopside, Hematite, Pyrite, Pyrrhotite, Rutile, Zircon, and most `.surf` variables) were driven to zero by the L1 penalty.

### Features that Increase Flotation Rate

| Feature | Coefficient | Hazard Ratio | 95% CI | p-value | Interpretation |
|---------|------------|--------------|--------|---------|----------------|
| **Apatite** | +0.55 | 1.74 | [0.53, 0.57] | <0.005 | The ore mineral. Higher apatite content strongly increases flotation rate. Each percentage point increase in apatite raises the hazard by 74%, meaning faster recovery. |
| **ECD** | +0.02 | 1.02 | [0.02, 0.03] | <0.005 | Equivalent Circle Diameter (particle size). Larger particles float slightly faster, consistent with flotation kinetics theory where larger particles have greater bubble attachment probability. |

### Features that Decrease Flotation Rate

| Feature | Coefficient | Hazard Ratio | 95% CI | p-value | Interpretation |
|---------|------------|--------------|--------|---------|----------------|
| **Phlogopite.surf** | −0.99 | 0.37 | [-1.08, -0.90] | <0.005 | Phlogopite (a phyllosilicate mica) surface exposure is the strongest inhibitor of flotation. Each unit increase in phlogopite surface fraction reduces the hazard by 63%. Mica surfaces are hydrophilic and create a wetting barrier that prevents bubble attachment. |
| **Quartz.surf** | −0.52 | 0.59 | [-0.62, -0.43] | <0.005 | Quartz surface exposure reduces flotation rate by 41%. Quartz is a classic non-floatable gangue mineral — its surface is naturally hydrophilic. |
| **Albite.surf** | −0.52 | 0.60 | [-0.69, -0.36] | <0.005 | Albite (feldspar) surface exposure reduces flotation rate by 40%, similar to quartz. Both are silicate gangue minerals. |
| **Orthoclase.surf** | −0.30 | 0.74 | [-0.48, -0.12] | <0.005 | Orthoclase (another feldspar) surface reduces flotation rate by 26%, a weaker effect than albite. |
| **Sanidine.surf** | −0.21 | 0.81 | [-0.43, -0.00] | <0.005 | Sanidine (high-temperature feldspar) surface has a modest inhibitory effect (19% reduction). |

## Key Physical Insights

### 1. Surface Mineralogy Dominates Over Bulk Composition

A critical finding: the model retained surface exposure fractions (`.surf` features) as the significant gangue predictors while eliminating corresponding bulk composition variables. This confirms that **what is on the particle surface determines flotation behavior**, not what is inside the particle.

### 2. Apatite (Ore) vs. Gangue Surface: A Clear Separation

Particles rich in apatite float rapidly (hazard ratio 1.74), while particles with exposed gangue minerals (phlogopite, quartz, albite, orthoclase) float slowly or not at all. This creates a natural separation basis for flotation circuit optimization.

### 3. Phlogopite is the Most Problematic Gangue Mineral

Phlogopite surface exposure has the strongest negative effect (coefficient −0.99). Mica minerals pose a particular challenge in phosphate flotation because:
- They are platy and tend to be entrained mechanically
- Their basal surfaces are hydrophilic, resisting collector adsorption
- When locked with apatite, they reduce composite particle floatability

### 4. Particle Size Effect

The positive ECD coefficient (+0.023) confirms that larger particles in this size fraction float slightly faster. While the effect is small per micron, it is highly significant (z = 61.5) across the particle size distribution.

## Methodology Notes

### Data Balancing
The original dataset (796,130 rows) was balanced to 300,000 rows with 60,000 samples per class (CA20, CB20, CC20, CD20, TD20). This prevents the model from being biased toward overrepresented classes.

### Feature Selection Pipeline
1. **Constant removal**: Dropped Clinochlore and Clinochlore.surf (zero variance)
2. **VIF filtering**: Removed 27 features with VIF > 10, retaining 24 features
3. **LASSO selection**: L1 regularization drove 18 feature coefficients to approximately zero

### Bootstrap Stability
All 1000 bootstrap iterations converged. The narrow confidence intervals (particularly for Apatite, Phlogopite.surf, and Quartz.surf) demonstrate that these features are robust predictors across different data subsets.
