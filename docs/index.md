# SurvBregDiv: Transfer Learning for Time-to-Event Modelling via Bregman Divergence

**SurvBregDiv** provides flexible and efficient tools for integrating
external risk scores into time-to-event models (including Cox
proportional hazards models and nested case–control models), while
properly accounting for population heterogeneity. The package enables
robust estimation, improved predictive performance, and user-friendly
workflows for modern survival analysis.

## Key Features

- **Integration of External Risk Scores**  
  Seamlessly incorporate predictions from *any* external risk model into
  your analysis pipeline.

- **Population Heterogeneity Adjustment**  
  Adjusts for distributional differences between the external model’s
  source population and the target study population.

- **Multiple Data Structures and Modelling Settings**  
  Supports both Cox proportional hazards models and nested case–control
  designs, with functionality for low-dimensional inference and
  high-dimensional feature selection.

- **Improved Estimation and Prediction**  
  Empirically demonstrates efficiency gains in parameter estimation and
  enhanced predictive accuracy.

- **Built-In Cross-Validation**  
  Automated procedures for tuning and penalization parameter selection.

## Installation

> **Note:** This package is currently under active development. Please
> report any issues or unexpected behavior.  
> Requires **R ≥ 4.0.0**.

Install from CRAN:

``` R
install.packages("SurvBregDiv")
```

Or install the development version from GitHub:

``` R
require("devtools")
require("remotes")
remotes::install_github("UM-KevinHe/SurvBregDiv")
```

## Detailed Tutorial

Full package documentation and parameter explanations:
[here](https://um-kevinhe.github.io/SurvBregDiv/)

## Getting Help

If you encounter problems or bugs, please contact us:

- <ybshao@umich.edu>
- <junyiqiu@umich.edu>
- <kevinhe@umich.edu>
