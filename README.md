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

Machine learning (ML) models have become increasingly popular in psychological research, often yielding remarkable results. However, they have also faced growing criticism in recent times. The discrepancy between proposed claims and reality is commonly attributed to the quality of psychological datasets, which tend to be small and subject to imprecise measurement. In this simulation study, we examined the data requirements necessary for complex ML algorithms to perform well. We compared the performance of Elastic Net Regressions and Gradient Boosting Machines (GBM) under various conditions: (a) sample size, (b) number of irrelevant predictors, (c) predictor reliability, (d) effect size, and (e) nature of the data generating model (i.e., linear vs. non-linear effects). We investigated whether the models achieved the maximum predictive performance (i.e., recovery of the simulated effects) and whether the respective models were able to accurately distinguish between relevant and irrelevant predictors. Our results showed that Elastic Net Regressions with pre-specified interaction terms outperformed GBM models in almost all conditions, but still only achieved the maximum possible performance under optimal conditions (\textit{N} = 1,000, perfectly reliable predictors, predominantly linear effects, and an exceptionally large effect size of \textit{R}² = .80), which are rarely met in psychological research. In general, we stress that data quality fundamentally limits the performance of ML models in a similar way to more traditional regression analyses.

---

## Folder Structure

project/\
  ├── utils/ # Utility functions for simulation and analysis\
  │ └── [utility scripts, e.g., anaylsisTools.R, fitENET.R, setParameters.R]\
  ├── supplement/ # code and plots for supplementary analyses, etc.\
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
  ├── 05b_plotMainPVI_EnetGBM.R\
  ├── 06a_analyseInterPVI_ENETw.R\
  ├── 06b_plotInterPVI_ENETw.R\
  ├── 06c_analyseInterStrength_GBM.R\
  ├── 06d_plotInterStrength_GBM.R\
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
     - Performs ANOVA to assess model performance and plots results of the generalized $\eta^2$. (`04a_analyseR2_ENETGBM.R`)
     - Prepares data for visualizations (creates `rSquaredData_stats.rda`) and generates an initial overview plot for $R^2$ (+ overfit). (`04b_prepareR2_EnetGBM.R`)
     - Creates result plots. (`04c_plotR2_EnetGBM.R`)
     
5. **05[a-c]_analyse[.]_[.].R**:
   - **Purpose**: Detects and analyzes linear effects.
   - **Process**:
     - Prepares data for visualizations (creates `relFrequencyMeasures.rda`) and performs ANOVA for sensitivity and specificity. (`05a_analyseMainPVI_EnetGBM.R`)
     - Creates result plots. (`05b_plotMainPVI_EnetGBM.R`)

6. **06[a-d]_[.].R**:
   - **Purpose**: Analyzes and plots interaction effects.
   - **Process**:
     - Prepares ENET data for visualizations (creates `interENETinterMeasures.rda`) and performs ANOVA for sensitivity and specificity. (`06a_analyseInterPVI_ENETw.R`)
     - Creates result plots for the ENETinter. (`06b_plotInterPVI_ENETw.R`)
     - Prepares GBM data for visualizations (creates `hStatsPlottingData.rda`) and performs ANOVA for the H-statistic. (`06c_analyseInterStrength_GBM.R`)
     - Creates result plots for the H-statistic of the GBM. (`06d_plotInterStrength_GBM.R`)

---

### Utility Files

- Located in the `utils/` folder.
- Files and functions to facilitate...
  - ... overviewing simulation conditions and parameter setup (`setParameters.R`).
  - ... simulation of the data and analysis of the results (`simTools.R`, `analysisTools.R`)
  - ... model fitting (`fitENET.R`, `fitGBM.R`, `saveENET.R`, `saveGBM.R`)

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
   - Use the {04, 05, 06} scripts to replicate the results.

---

## Contact Information

If you have any questions or encounter issues, please contact:
- **Name**: Kim-Laura Speck
- **Email**: kim.speck@uni-kassel.de
- **Affiliation**: Universität Kassel

---