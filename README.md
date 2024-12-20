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

Machine Learning (ML) models have brought valuable methodological innovations to psychological research. However, studies applying non-parametric ML models to enhance the prediction of outcomes in a psychological context have yielded mixed results, with performance frequently comparable to that of traditional regression models. This discrepancy between aspiration and reality is commonly attributed to the nature and quality of psychological datasets, which tend to be small and subject to imprecise measurement in comparison to the larger, more robust datasets used in computational sciences, where ML models were developed. In this study, we examined the data requirements necessary for complex algorithms to perform well. In an extensive simulation study, we compared the performance of Elastic Net Regressions and Gradient Boosting Machines under systematically varied conditions: (a) sample size, (b) number of irrelevant predictors, (c) amount of measurement error, (d) effect sizes (i.e., proportion of explained variance), and (e) the nature of the data generating model (i.e., linear vs.  non-linear effects). We investigated two primary criteria: First, whether the models achieved the maximum predictive performance (i.e., recovery of the simulated systematic variance). Second, whether the respective models were able to accurately identify the correct predictors as important variables. The overall goal of this simulation study is to guide expectations regarding the achievable performance of sophisticated ML algorithms in psychological research and to identify scenarios where their application may obscure, rather than reveal, underlying relationships.

---

## Folder Structure

project/\
  ├── utils/ # Utility functions for simulation and analysis\
  │ └── [utility scripts, e.g., anaylsisTools.R, fitENET.R, setParameters.R]\
  ├── supplement/ # code for supplementary analyses, etc.\
  ├── results/ # Folder for fitted models and dependent measure files\
  │ ├── finalResults/ # Contains model results in .rds files\
  │   ├── dependentMeasures/ # Subfiles for dependent measures\
  ├── plots/ # Folder for result plots\
  │ ├── detectMains/ # Contains plots for the detection of linear effects\
  │ ├── detectInteractions/ # Contains plots for the detection of interaction effects\
  │ ├── ANOVAresults/ # Contains plots for the ANOVA results\
  ├── log/ # log files from data simulation and fitting\
  ├── info/ # Information about results file structures\
  │ └── resultVariablesOverview.Rmd # Metadata documentation\
  ├── data/ # Folder containing simulated data\
  │   └── [simulated data files; this data is not provided due to memory limitations]\
  ├── checkSimulation/ # Folder containing files for initial sanity checks of the simulation\
  ├── 01_simulateData.R # Script to simulate raw data\
  ├── 02_fitData.R # Script to fit models and save results\ 
  ├── 03_joinData.R # Script to merge data across conditions\
  ├── 04a_analyseR2_EnetGBM.R # Analyzes R² results\
  ├── 04b_prepareR2_EnetGBM.R # \
  ├── 04c_plotR2_EnetGBM.R # \
  ├── 05a_analyseMainPVI_EnetGBM.R # analyzes linear effects based on the PVI values\
  ├── 05b_analyseMainBeta_ENETw.R # \
  ├── 06a_analyseInterPVI_ENETw.R # \
  ├── 06b_analyseInterStrength_GBM.R # Detects and analyzes interaction effects\
  ├── MLsim.Rproj # R project\
  └── README.md # Documentation (this file)\
---

## Getting Started

1. **Dependencies**:
    - R (>= 4.3)
    - Required R packages: `mvtnorm`, `truncnorm`, `parallel`, `glmnet`, `gbm`, `caret`, `iml`, etc.

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
    4. Analyse data: {04, 05, 06} files

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
     - Outputs merged data files stored in the `results/finalResults/dependentMeasures` folder.

4. **04[a-c]_[.]_ENETGBM.R** and following files:
   - **Purpose**: Analyzes \(R^2\) results to determine the effect of experimental manipulations.
   - **Process**:
     - Performs ANOVA to assess model performance. (`04a_analyseR2_ENETGBM.R`)
     - Prepares data for visualizations. (`04b_prepareR2_EnetGBM.R`)
     - Creates result plors. (`04c_plotR2_EnetGBM.R`)
     
5. **05[a-c]_analyse[.]_[.].R**:
   - **Purpose**: Detects and analyzes linear effects.
   - **Process**:
     - Runs ANOVA for each condition and dependent variable.
     - Generates data for plots.

6. **analyze_interaction_effects.R**:
   - **Purpose**: Detects and analyzes interaction effects.
   - **Process**:
     - Similar to `analyze_linear_effects.R`, but focuses on interaction terms.
     - Creates visual summaries of results.

---

### Utility Files

- Located in the `utils/` folder.
- Functions to facilitate:
  - Simulation and parameter setup.
  - Efficient storage and retrieval of simulation parameters.

---

### Data Folder

- The `data/` folder contains raw simulated data produced by `simulate_data.R`.

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
1. Simulate data using `simulate_data.R`.
2. Fit models with `fit_data.R`.
3. Merge and analyze data:
   - Run `join_data.R`.
   - Use the `analyze_*` scripts to explore specific aspects of the results.

---

## Contact Information

If you have any questions or encounter issues, please contact:
- **Name**: [Your Name]
- **Email**: [Your Email Address]
- **Affiliation**: [Your Institution]

---