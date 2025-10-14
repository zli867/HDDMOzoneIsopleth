# Ozone Isopleth Model Suite

A collection of MATLAB tools for developing and analyzing ozone response surfaces using quadratic and logarithmic quadratic models with NOx and VOC precursor sensitivity analysis.

## Quadratic Models

### Full Quadratic Model
**File:** `QuadraticModel.m`

Develops a complete second-order polynomial model relating ozone concentrations to NOx and VOC emissions:

```
O3 = a0 + a1*NOx + a2*VOC + a3*NOx² + a4*NOx*VOC + a5*VOC²
```

### Reduced Quadratic Model
**File:** `QuadraticReduce.m`

A simplified variant that excludes the VOC quadratic term, accounting for the lower sensitivity of ozone formation to VOC in certain regimes:

```
O3 = a0 + a1*NOx + a2*VOC + a3*NOx² + a4*NOx*VOC
```

### Weighted Quadratic Model
**File:** `WeightedQuadraticModel.m`

An extension of the full quadratic model that applies differential weighting to data points based on sensitivity. Enables generation of weighted response surfaces and isopleths that emphasize high-sensitivity regions.

---

## Logarithmic Quadratic Models

### Full Log-Quadratic Model
**File:** `LogQuadraticModel.m`

Develops a second-order polynomial model in log-space, capturing nonlinear responses more effectively:

```
log(O3) = a0 + a1*NOx + a2*VOC + a3*NOx² + a4*NOx*VOC + a5*VOC²
```

### Reduced Log-Quadratic Model
**File:** `LogQuadraticReduce.m`

Simplified log-space model excluding the VOC quadratic term:

```
log(O3) = a0 + a1*NOx + a2*VOC + a3*NOx² + a4*NOx*VOC
```

### Weighted Log-Quadratic Model
**File:** `WeightedLogQuadratic.m`

Extends the full log-quadratic model with differential weighting for sensitivity-based isopleth generation.

### Alternative Log-Base Formulation
**File:** `LogQuadraticDiffBase.m`

Develops log-quadratic models using alternative logarithmic bases. While mathematically equivalent to the standard natural log approach, this variant provides flexibility for different analytical frameworks.

---

## Uncertainty Analysis

### Leave-One-Out Validation

**Quadratic Model:** `QuadraticWithholding.m`
**Log-Quadratic Model:** `LogQuadraticWithholding.m`

Both programs implement a data withholding (cross-validation) approach:
- Sequentially omits each data point
- Regenerates the ozone isopleth without the withheld data
- Calculates mean bias and standard deviation of prediction errors
- Provides quantitative assessment of model robustness and uncertainty

---

## Trajectory and Sensitivity Analysis

### Sensitivity Trajectory Analysis

**File:** `SensitivityTrajectory.m`

Performs comprehensive sensitivity analysis combining observational data and model predictions:
- Estimates empirical ozone sensitivities from observations
- Compares quadratic and log-quadratic model sensitivities
- Generates four-panel comparison plots showing:
  - Modeled vs. observed O3
  - Model sensitivity to NOx (dO3/dNOx)
  - Model sensitivity to VOC (dO3/dVOC)
  - Cross-model comparison metrics

### Comparative Sensitivity Trajectory Analysis

**File:** `SensitivityTrajectoryCompare.m`

Streamlined two-panel visualization overlaying results from both model types:
- Quadratic and log-quadratic results displayed together
- Direct comparison of model performance and sensitivities
- Compact format for rapid inter-model evaluation

### Regression Line Analysis

**File:** `RegressionLine.m`

Scatter plot analysis of least-squares model performance:
- Plots CMAQ-DDM input precursor data against isopleth model outputs
- Calculates regression statistics and residuals
- Assesses model fit quality and systematic biases

---

## Model Comparison and Validation

### HDDM vs. Observation-Based Models

**File:** `compareQuadraticEmp.m`

Comparative analysis between two independent modeling approaches:
- Contrasts quadratic-form HDDM (Hierarchical Deterministic Dynamic Model) results with observation-derived models
- Evaluates differences in ozone isopleth construction methodologies
- Identifies strengths and limitations of each approach

### CMAQ-HDDM vs. Fitted Quadratic/Log-Quadratic Models

**File:** `TrajcompareCMAQ.m`

Validates fitted model performance against full CMAQ-HDDM simulations:
- Compares predictions from quadratic or log-quadratic models trained on CMAQ-HDDM data
- Assesses model accuracy in reproducing original CMAQ results

### Observation-Derived vs. CMAQ-HDDM Sensitivities

**File:** `TrajCompare.m`

Cross-validates three data sources for ozone isopleth sensitivity:
- Compares ozone observations with CMAQ-HDDM-derived ozone isopleth predictions
- Assesses consistency of dO3/dNOx and dO3/dVOC across different methodologies

---
