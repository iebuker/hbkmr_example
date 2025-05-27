# Bayesian Kernel Machine Regression for Heteroscedastic Health Outcome Data

This repository accompanies the paper:

**_Bayesian Kernel Machine Regression for Heteroscedastic Health Outcome Data_**

In this work, we introduce a Bayesian framework for modeling heteroscedastic health outcome data using kernel machine regression. The method is implemented in [NIMBLE](https://r-nimble.org) and includes tools for model fitting, diagnostics, and visualizing both residual variation and the estimated $h(\cdot)$ function.

## Getting Started

To explore the method using a simulated dataset, open:

`overview.html`

This document provides a step-by-step demonstration of:

- Simulation of heteroscedastic data
- Model fitting using NIMBLE
- Visualization of Bayesian residuals to inform the variance model 
- Posterior inference using helper functions 
- Visualization of results

## Key Files

- `SimHDData.R`: Simulates heteroscedastic health outcome data
- `makeK.R`: Constructs the kernel matrix
- `checkResidualVar.R`, `helper_checkResidualVar.R`: Functions for examining heteroscedasticity in the residuals
- `createPlots.R`, `helper_createPlots.R`: Functions for plotting the estimated \( h(\cdot) \) function
- `overview.Rmd`: R Markdown source for the overview
- `overview.html`: Rendered HTML walkthrough of the example
