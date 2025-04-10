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

---

### 2. Choose a Simulation Scenario

Go to the relevant folder under `simulation/`:

- `simResults_scenario`: frequent exposure and outcome
- `simResults_scenarioER`: rare exposure
- `simResults_scenarioOR`: rare outcome

---

### 3. Run Simulations for Each Method

**Step 1: Navigate to the Method Folder**

- Inside each scenario folder, go to the relevant method folder.

**Step 2: Run the Script**

- Open and run the script to generate results.
- To control the number of simulation iterations, set:
  - `a`
- Total iterations = `a`, i.e., if `a` is set to 10, the number of simulation iterations is 10.

**Step 3: View the Results**

- Individual results are saved to `result/results_METHOD.X.Rds` files in the method folder. `METHOD` = the method folder name; `X` = the simulation iteration index.
  - Some individual results contain `NULL` values due to the specific simulation datasets being incompatible with the applied method.
- These are aggregated into the `RD_METHOD` object in memory. `METHOD` = the method folder name.
- To save the combined output, manually uncomment the final two lines of the code.

---

## 📝 Notes

- The output `RD_METHOD` is created in memory and not saved unless modified. (`METHOD` = the method folder name.)
- Keep the file structure unchanged unless necessary.
  
---

## 📄 License

This project is licensed under the GPL-3.0 license.
