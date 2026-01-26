# Plot Cross-Validation Results vs Eta (Biometrics-style)

Plots cross-validation performance across eta values for `cv.coxkl`,
`cv.coxkl_ridge`, or `cv.coxkl_enet` objects in the Biometrics figure
style. It displays the CV performance curve, a baseline reference at
eta=0, and marks the optimal eta.

## Usage

``` r
cv.plot(object, line_color = "#7570B3", baseline_color = "#1B9E77", ...)
```

## Arguments

- object:

  A fitted cross-validation result of class `"cv.coxkl"`,
  `"cv.coxkl_ridge"`, or `"cv.coxkl_enet"`.

- line_color:

  Color for the CV performance curve. Default is `"#7570B3"`.

- baseline_color:

  Color for the external baseline line. Default is `"#1B9E77"`.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot` object showing cross-validation performance versus `eta`.

## Examples

``` r
data(ExampleData_lowdim)
train_dat_lowdim <- ExampleData_lowdim$train
test_dat_lowdim <- ExampleData_lowdim$test
beta_external_lowdim <- ExampleData_lowdim$beta_external_good

eta_list <- generate_eta(method = "exponential", n = 50, max_eta = 30)
cv.coxkl_est <- cv.coxkl(z = train_dat_lowdim$z,
                         delta = train_dat_lowdim$status,
                         time = train_dat_lowdim$time,
                         beta = beta_external_lowdim,
                         etas = eta_list,
                         criteria = "V&VH",
                         seed = 1)
#> Warning: Stratum not provided. Treating all data as one stratum.
cv.plot(cv.coxkl_est)

```
