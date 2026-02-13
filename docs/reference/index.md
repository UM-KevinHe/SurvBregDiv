# Package index

## Individual-level Internal–External Integration

- [`cox_indi()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cox_indi.md)
  : Cox Proportional Hazards Model Integrated with External
  Individual-level Information
- [`cv.cox_indi()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.cox_indi.md)
  : Cross-Validated cox_indi to Tune etas

## Low-Dimensional Bregman Divergence based Integrated Cox Model

- [`coxkl()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl.md)
  : Cox Proportional Hazards Model with KL Divergence for Data
  Integration
- [`cv.coxkl()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.coxkl.md)
  : Cross-Validated Cox–KL to Tune the Integration Parameter (eta)
- [`coxkl_ties()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_ties.md)
  : Cox Proportional Hazards Model with KL Divergence for Data
  Integration (Ties Handling)
- [`cv.coxkl_ties()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.coxkl_ties.md)
  : Cross-Validated Cox–KL with Ties Handling to Tune the Integration
  Parameter (eta)
- [`cox_MDTL()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cox_MDTL.md)
  : Cox Proportional Hazards Model with Mahalanobis Distance Transfer
  Learning
- [`cv.cox_MDTL()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.cox_MDTL.md)
  : Cross-Validation for Cox MDTL Model

## High-Dimensional Bregman Divergence based Integrated Cox Model

- [`coxkl_ridge()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_ridge.md)
  : Cox Proportional Hazards Model with Ridge Penalty and External
  Information
- [`coxkl_enet()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_enet.md)
  : Cox Proportional Hazards Model with KL Divergence and Elastic Net
  Penalty
- [`cv.coxkl_ridge()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.coxkl_ridge.md)
  : Cross-Validation for CoxKL Ridge Model (Tuning Eta and Lambda)
- [`cv.coxkl_enet()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.coxkl_enet.md)
  : Cross-Validation for CoxKL Model with Elastic Net & Lasso Penalty
- [`cox_MDTL_ridge()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cox_MDTL_ridge.md)
  : Cox MDTL with Ridge Regularization
- [`cox_MDTL_enet()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cox_MDTL_enet.md)
  : Fit Cox Model with Multi-Domain Transfer Learning and Elastic Net
  Penalty
- [`cv.cox_MDTL_ridge()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.cox_MDTL_ridge.md)
  : Cross-Validation for Cox MDTL with Ridge Regularization
- [`cv.cox_MDTL_enet()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.cox_MDTL_enet.md)
  : Cross-Validation for Cox MDTL with Elastic Net Regularization

## Variable Importance/Stability Selection for High-Dimensional Models

- [`variable_importance()`](https://um-kevinhe.github.io/SurvBregDiv/reference/variable_importance.md)
  : Bootstrap Variable Importance via Selection Frequency
- [`coxkl_enet.StabSelect()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_enet.StabSelect.md)
  : Stability Selection for KL-Integrated Cox Elastic-Net Models
- [`cox_MDTL_enet.StabSelect()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cox_MDTL_enet.StabSelect.md)
  : Stability Selection for MDTL-Integrated Cox Elastic-Net Models

## Bagging for High-Dimensional Models

- [`coxkl_enet_bagging()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_enet_bagging.md)
  : Bagging for KL-Integrated Cox Elastic-Net Models
- [`cox_MDTL_enet_bagging()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cox_MDTL_enet_bagging.md)
  : Bagging for MDTL-Integrated Cox Elastic-Net Models

## Low-Dimensional KL Divergence based Integrated Conditional Logistic Model

- [`clogitkl()`](https://um-kevinhe.github.io/SurvBregDiv/reference/clogitkl.md)
  : Conditional Logistic Regression with KL Divergence (CLR-KL)
- [`cv.clogitkl()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.clogitkl.md)
  : Cross-Validated Conditional Logistic Regression with KL Integration

## High-Dimensional KL Divergence based Integrated Conditional Logistic Model

- [`clogitkl_enet()`](https://um-kevinhe.github.io/SurvBregDiv/reference/clogitkl_enet.md)
  : Conditional Logistic Regression with KL Divergence and Elastic Net
  Penalty (CLR-KL-ENet)
- [`cv.clogitkl_enet()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.clogitkl_enet.md)
  : Cross-Validated CLR-KL with Elastic Net Penalty

## Multi-Source Integration

- [`coxkl_enet.multi()`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_enet.multi.md)
  : Multi-Source Integration for KL-Integrated Cox Elastic-Net Models

## Plotting Functions

- [`cv.plot()`](https://um-kevinhe.github.io/SurvBregDiv/reference/cv.plot.md)
  : Plot Cross-Validation Results vs Eta
- [`plot(`*`<coxkl>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.coxkl.md)
  : Plot Validation Results for coxkl Object
- [`plot(`*`<coxkl_ridge>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.coxkl_ridge.md)
  : Plot Validation Results for coxkl_ridge Object
- [`plot(`*`<coxkl_enet>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.coxkl_enet.md)
  : Plot Validation Results for coxkl_enet Object
- [`plot(`*`<cox_MDTL>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.cox_MDTL.md)
  : Plot Validation Results for Cox_MDTL Object
- [`plot(`*`<cox_MDTL_enet>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.cox_MDTL_enet.md)
  : Plot Validation Results for cox_MDTL_enet Object
- [`plot(`*`<cox_MDTL_ridge>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.cox_MDTL_ridge.md)
  : Plot Validation Results for cox_MDTL_ridge Object
- [`plot(`*`<StabSelect>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.StabSelect.md)
  : Plot Stability Selection Path
- [`plot(`*`<variable_importance>`*`)`](https://um-kevinhe.github.io/SurvBregDiv/reference/plot.variable_importance.md)
  : Plot Variable Importance (Selection Frequency)

## Utilities

- [`generate_eta()`](https://um-kevinhe.github.io/SurvBregDiv/reference/generate_eta.md)
  : Generate a Sequence of Tuning Parameters (eta)

## Datasets

- [`ExampleData_lowdim`](https://um-kevinhe.github.io/SurvBregDiv/reference/ExampleData_lowdim.md)
  : Example low-dimensional survival data
- [`ExampleData_highdim`](https://um-kevinhe.github.io/SurvBregDiv/reference/ExampleData_highdim.md)
  : Example high-dimensional survival data
- [`ExampleData_cc`](https://um-kevinhe.github.io/SurvBregDiv/reference/ExampleData_cc.md)
  : Example Data for Conditional Logistic Regression
- [`ExampleData_cc_highdim`](https://um-kevinhe.github.io/SurvBregDiv/reference/ExampleData_cc_highdim.md)
  : Example high-dimensional matched case-control data
- [`ExampleData_indi`](https://um-kevinhe.github.io/SurvBregDiv/reference/ExampleData_indi.md)
  : Example internal/external Cox individual-level data
