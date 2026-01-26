# Stability Selection for Cox-KL Enet Model

Stability Selection for Cox-KL Enet Model

## Usage

``` r
coxkl_enet.StabSelect(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  etas,
  alpha = 1,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = 0.1,
  nfolds = 5,
  cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
  c_index_stratum = NULL,
  message = FALSE,
  seed = NULL,
  B = 50,
  fraction_sample = 0.5,
  ...
)
```

## Arguments

- z:

  Numeric matrix of covariates. Rows represent individuals and columns
  represent predictors.

- delta:

  Numeric vector of event indicators (1 = event, 0 = censored).

- time:

  Numeric vector of observed times (event or censoring).

- stratum:

  Optional numeric or factor vector indicating strata. If `NULL`, all
  subjects are assumed to be in the same stratum.

- RS:

  Optional numeric vector or matrix of external risk scores. If not
  provided, `beta` must be supplied.

- beta:

  Optional numeric vector of external coefficients (length equal to
  `ncol(z)`). If provided, it is used to compute external risk scores.
  If not provided, `RS` must be supplied.

- etas:

  Numeric vector of candidate `eta` values to be evaluated.

- alpha:

  Elastic-net mixing parameter in \\(0,1\]\\. Default is `1` (lasso
  penalty).

- lambda:

  Optional numeric vector of lambda values. If `NULL`, a path is
  generated automatically.

- nlambda:

  Integer. Number of lambda values to generate if `lambda` is NULL.
  Default is 100.

- lambda.min.ratio:

  Numeric. Ratio of min/max lambda. Default depends on sample size vs
  dimension (0.05 if n \< p, else 1e-03).

- nfolds:

  Integer. Number of cross-validation folds. Default is `5`.

- cv.criteria:

  Character string specifying the cross-validation criterion for
  selecting both `eta` and `lambda`. Choices are:

  - `"V&VH"` (default): V&VH loss.

  - `"LinPred"`: Loss based on cross-validated linear predictors.

  - `"CIndex_pooled"`: Pooled C-Index.

  - `"CIndex_foldaverage"`: Average C-Index across folds.

- c_index_stratum:

  Optional stratum vector. Required only when `cv.criteria` is set to
  `"CIndex_pooled"` or `"CIndex_foldaverage"`, and a stratified C-index
  needs to be computed while the fitted model is non-stratified. Default
  is `NULL`.

- message:

  Logical. Whether to print progress messages. Default is `FALSE`.

- seed:

  Optional integer. Random seed for reproducible fold assignment.

- B:

  Integer. Number of bootstrap/subsampling replicates. Default is 50.

- fraction_sample:

  Numeric. Fraction of data to use for subsampling in each replicate.
  Default is 0.5.

- ...:

  Additional arguments passed to the underlying fitting function
  [`coxkl_enet`](https://um-kevinhe.github.io/SurvBregDiv/reference/coxkl_enet.md).

## Value

An object of class `StabSelect`, which is a list containing:

- stability_path:

  A matrix (n_vars x n_lambda) of selection probabilities.

- lambda:

  The global lambda sequence used.

## Examples

``` r
if (FALSE) { # \dontrun{
data(ExampleData_highdim)

# Plot with different thresholds without re-running
plot(coxkl.StabSelect, threshold = 0.6)
plot(coxkl.StabSelect, threshold = 0.8)
} # }
```
