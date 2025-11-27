# MLsim

This repository contains the code and resources for a simulation study designed to investigate the performance of ML models on psychological data. The code is structured to enable reproducibility and transparency of the results.

---

## Table of Contents
1. [Overview](#overview)
2. [Folder Structure](#folder-structure)
3. [Getting Started](#getting-started)
4. [File Descriptions](#file-descriptions)
5. [Usage](#usage)
6. [Contact Information](#contact-information)

---

## Overview

Although Machine Learning (ML) methods are gaining popularity in psychological research, the debate about their usefulness ranges from hype to disillusionment. The discrepancy between the hopes placed in ML methods and the empirical reality is often attributed to the quality of psychological datasets, which tend to be small and subject to imprecise measurement. In this simulation study, we examined the data requirements necessary for ML methods to perform well. We compared the performance of Elastic Net Regressions with and without prespecified interactions, Random Forests and Gradient Boosting Machines for different data-generating processes (including either interaction, stepwise, or piecewise linear effects) and under various conditions: (a) sample size, (b) number of irrelevant predictors, (c) predictor reliability, (d) effect size, and (e) nature of the data generating model (i.e., linear vs. non-linear effects). We investigated whether the models achieved the highest level of predictive performance attainable under the given simulated conditions. There were two main takeaways of our results: First, the maximum possible predictive performance was only achieved under optimal simulation conditions (N = 1,000, perfectly reliable predictors, predominantly linear effects, and an exceptionally large effect size of R² = .80), which are arguably rarely met in psychological research. Second, each ML model outperformed the others under certain conditions, but none was consistently superior or entirely robust to suboptimal data characteristics. We stress that data quality fundamentally limits predictive performance and discuss the interpretation of comparisons between flexible ML models and simpler (regularized linear) baselines in psychological research.

---

## Folder Structure

project/\
  ├── 01_simulateData.R # Script to simulate raw data\
  ├── 02_fitData.R # Script to fit models and save results\
  ├── 03_joinData.R # Script to merge data across conditions\
  ├── 04_analyseR2_ANOVA.R # Analyzes R² results\
  ├── 05_plotR2_plotOverfit.R # \
  ├── MLsim.Rproj # R project\
  ├── README.md # Documentation (this file)\
  ├── utils/ # Utility functions for simulation and analysis\
  │ └── [utility scripts, e.g., anaylsisTools.R, fitENET.R, setParameters.R]\
  ├── onlineMaterial/ # code and plots for supplementary analyses, etc.\
  ├── results/ # Folder for fitted models and dependent measure files\
  │ ├── dependentMeasures/
  │ ├── pwlinear/dependentMeasures/ # Contains model results for respective DGP\
  │ ├── nonlinear3/dependentMeasures/ # Contains model results for respective DGP\
  │ ├── inter/dependentMeasures/ # Contains model results for respective DGP\
  ├── plots/ # Folder for result plots\
  │ ├── ANOVAresults/ # Contains plots for the ANOVA results\
  │ ├── hyperParamter/ # Contains histograms for hyperparameter choices\
  │ ├── ... paper png-files # Contains plots that made it into the paper\
  ├── info/ # Information about results file structures\
  │ └── resultVariablesOverview.Rmd # Metadata documentation\
  ├── data/ # Folder containing simulated data\
  │   └── [simulated data files; this data is not provided due to memory limitations]\
  └── log/ # log files from data simulation and fitting\

---

## Getting Started

1. **Dependencies**:
    - R (>= 4.5)
    - Required R packages: `mvtnorm`, `truncnorm`, `parallel`, `glmnet`, `gbm`, `caret`, etc.

2. **Installation**:
    Clone this repository:
    ```bash
    git clone https://github.com/kimSpeck/MLsim.git
    ```

3. **Usage**:
    Run the scripts in the specified order:
    1. Simulate data: 01_simulateData.R
    2. Fit models to the simulated data: 02_fitData.R
    3. Merge data: 03_joinData.R
    4. Analyse data: {04, 05, ...} files

---

## File Descriptions

### Main Scripts

1. **01_simulateData.R**:
   - **Purpose**: Simulates raw data using a regression model, sampling from a multivariate normal distribution.
   - **Key Features**:
     - uses parameter specifications from `utils/setParamaters.R`
     - Outputs raw data stored in the `data/` folder.

2. **02_fitData.R**:
   - **Purpose**: Fits models to the simulated data and stores results for each simulated condition.
   - **Output**: 
     - uses utility functions from the `utils/` folder (e.g., `fitENET.R`, `saveENET.R`, ...) 
     - Results saved as `.rds` files in the `results/` folder.

3. **03_joinData.R**:
   - **Purpose**: Merges data from all simulated conditions for analysis.
   - **Features**:
     - Creates subfiles for dependent measures to avoid RAM overflow when analysing results.
     - Outputs merged data files stored in the respective DGP related results folder (e.g., `results/inter/dependentMeasures`) and in the `results/dependentMeasures` folder.

4. **04_analyseR2_ANOVA.R**:
   - **Purpose**: Analyzes \(R^2\) results to determine the effect of experimental manipulations.
   - **Process**:
     - Performs ANOVA to assess model performance and plots results of the generalized $\eta^2$. 
     - Prepares data for visualizations and generates ANOVA results plots. 
     
5. **05_plotR2_plotOverfit.R**:
    - **Purpose**: Plot \(R^2\) and overfit result graphics
    - **Process**:
      - Plots results graphics for \(R^2_{test}\) and Overfit (paper-style).
      - Plots results graphics for the full simulation design (provided in online material).

---

### Utility Files

- Located in the `utils/` folder.
- Files and functions to facilitate...
  - ... overviewing simulation conditions and parameter setup (`setParameters.R`).
  - ... simulation of the data and analysis of the results (`simTools.R`, `analysisTools.R`)
  - ... model fitting (`fitENET.R`, `fitGBM.R`, `fitRF.R`, `saveENET.R`, `saveGBM.R`, `saveRF.R`)

---

### Data Folder

- The `data/` folder contains raw simulated data produced by `01_simulateData.R`.

---

### Results Folder

- The `results/` folder stores:
  - Fitted model outputs.
  - Subfiles for dependent measures, categorized for efficient analysis.

---

### Info Folder

- The `info/` folder contains:
  - Metadata regarding the structure of results files.
  - Documentation on stored variables and formats.

---

## Usage

Follow the script order for reproducibility:
1. Simulate data using `01_simulateData.R`.
2. Fit models with `02_fitData.R`.
3. Merge and analyze data:
   - Run `03_joinData.R`.
   - Use the {04, 05, ...} scripts to replicate the results.

---

## Contact Information

If you have any questions or encounter issues, please contact:
- **Name**: Kim-Laura Speck
- **Email**: kim.speck@uni-kassel.de
- **Affiliation**: Universität Kassel

---