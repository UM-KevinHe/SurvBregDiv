# Bagging for cv.cox_MDTL_enet

Implements bootstrap aggregation (bagging) for the `cv.cox_MDTL_enet`
model. It generates `B` bootstrap samples from the internal data, fits
the Cross-Validated Cox MDTL Elastic Net model on each sample, and
aggregates the resulting coefficients.

## Usage

``` r
cox_MDTL_enet_bagging(
  z,
  delta,
  time,
  stratum = NULL,
  beta = NULL,
  vcov = NULL,
  etas,
  alpha = 1,
  B = 100,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = ifelse(nrow(z) < ncol(z), 0.01, 1e-04),
  nfolds = 5,
  cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
  c_index_stratum = NULL,
  message = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- z:

  Matrix of predictors (n x p).

- delta:

  Vector of event indicators.

- time:

  Vector of survival times.

- stratum:

  Vector indicating the stratum.

- beta:

  Vector of fixed external coefficients (length p). This is prior
  information and is **not** resampled.

- vcov:

  Optional weighting matrix (p x p).

- etas:

  Sequence of eta values (transfer learning weights) to tune.

- alpha:

  The Elastic Net mixing parameter, with \\0 \le \alpha \le 1\\.
  `alpha=1` is the lasso penalty, and `alpha=0` the ridge penalty.
  Default is 1.0.

- B:

  Integer. Number of bootstrap replicates. Default is 100.

- lambda:

  Optional user-supplied lambda sequence.

- nlambda:

  Number of lambda values.

- lambda.min.ratio:

  Ratio of min/max lambda.

- nfolds:

  Number of CV folds for the inner cross-validation.

- cv.criteria:

  Cross-validation criteria.

- c_index_stratum:

  Stratum vector for C-index calculation (if different from model
  stratum).

- message:

  Logical. If TRUE, shows a progress bar.

- seed:

  Integer. Seed for reproducibility.

- ...:

  Additional arguments passed to `cv.cox_MDTL_enet`.

## Value

An object of class `"cox_MDTL_bagging"` containing:

- `best_beta`: The averaged coefficient vector across all valid
  bootstrap replicates.

- `all_betas`: A matrix (p x B) of coefficients from each bootstrap
  replicate.

- `B`: Number of requested replicates.

- `valid_replicates`: Number of replicates that successfully converged.

- `seed`: The seed used.

## Examples

``` r
# \donttest{
data(ExampleData_highdim)
train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external

etas <- generate_eta(method = "exponential", n = 10, max_eta = 10)
bagging_res <- cox_MDTL_enet_bagging(
  z = train_dat_highdim$z,
  delta = train_dat_highdim$status,
  time = train_dat_highdim$time,
  stratum = train_dat_highdim$stratum,
  beta = beta_external_highdim,
  vcov = NULL,
  etas = etas,
  alpha = 0.5, # Elastic Net mixing
  B = 5,
  cv.criteria = "CIndex_pooled",
  message = TRUE,
  seed = 123
)
#> Starting Bagging (B = 5 ) for cv.cox_MDTL_enet:
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
# }
```
