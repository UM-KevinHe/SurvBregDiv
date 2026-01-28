# Plot Cross-Validation Results vs Eta

Plots cross-validation performance across eta values for `cv.coxkl`,
`cv.coxkl_ridge`, `cv.coxkl_enet`, `cv.cox_MDTL`, `cv.cox_MDTL_ridge`,
`cv.cox_MDTL_enet`, `cv.clogitkl`, or `cv.clogitkl_enet` objects in a
Biometrics-style figure. It displays the cross-validated performance
curve, a baseline reference at `eta = 0`, and marks the optimal `eta`.

## Usage

``` r
cv.plot(object, line_color = "#7570B3", baseline_color = "#1B9E77", ...)
```

## Arguments

- object:

  A fitted cross-validation result object.

- line_color:

  Color for the CV performance curve. Default is `"#7570B3"`.

- baseline_color:

  Color for the baseline line. Default is `"#1B9E77"`.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot` object.
