# hdPS-proxy-select-codes
Code for replicating the simulation study and real-world analysis. Includes all data and scripts for model implementation, performance evaluation, and visualization.

---

## 📁 Folder Structure

```
real_data_analysis/
├── data/                         # Real-World Data
├── RealDataAnalysis.Rmd          # Real-Data-Analysis code
simulation/
├── simData/                      # Simulated datasets
├── simResults_scenario/          # Frequent exposure & outcome scenario
├── simResults_scenarioER/        # Rare exposure scenario
├── simResults_scenarioOR/        # Rare outcome scenario
```

### 📂 Mapping of Folder Names

Each folder of name prefixed with `simResults_` contains subfolders for methods. The folder names correspond to the following methods used in the study. Each folder contains code for a specific method.

| **Folder Name**        | **Method**                                                |
|-------------------------|------------------------------------------------------------|
| `1_kitchen`            | Kitchen sink model                                         |
| `2_bross`              | hdPS using Bross formula                                   |
| `3_hybrid`             | Hybrid of hdPS and LASSO                                   |
| `4_lasso`              | Least Absolute Shrinkage and Selection Operator (LASSO)    |
| `5_enet`               | Elastic Net                                                |
| `6_rf`                 | Random Forest                                              |
| `7_xgboost`            | XGBoost                                                    |
| `8_forward`            | Stepwise - Forward Elimination                             |
| `9_backward`           | Stepwise - Backward Elimination                            |
| `10_ga`                | Genetic Algorithm                                          |


---

## ⚙️ How to Run Simulations

### 1. Setup

- **R required** (version ≥ 4.1 recommended)
- Install dependencies:
  ```r
  install.packages(c(
  "autoCovariateSelection", # hdPS analysis
  "glmnet",                 # LASSO and elastic-net
  "randomForest",           # Random forest algorithm for classification/regression  
  "xgboost",                # Extreme Gradient Boosting
  "leaps"                   # Best subset selection for regression models  
  "GA",                     # Genetic algorithm optimization
  "MASS",                   # Modern applied statistics
  "WeightIt",               # PS weighting
  "cobalt",                 # Balance tables and plots
  "ggplot2",                # Data visualization and plotting
  "dplyr"                   # Data manipulation
  ))
  ```

