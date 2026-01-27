# Bagging for cv.coxkl_enet

Bagging for cv.coxkl_enet

## Usage

``` r
coxkl_enet_bagging(
  z,
  delta,
  time,
  stratum = NULL,
  RS = NULL,
  beta = NULL,
  etas,
  alpha = 1,
  B = 100,
  lambda = NULL,
  nlambda = 100,
  lambda.min.ratio = ifelse(nrow(z) < ncol(z), 0.05, 0.001),
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

- RS:

  Vector of external risk scores. Resampled during bagging.

- beta:

  Vector of fixed external coefficients. Not resampled.

- etas:

  Sequence of eta values.

- alpha:

  Elastic net mixing parameter.

- B:

  Number of bootstrap replicates.

- lambda:

  Optional lambda sequence.

- nlambda:

  Number of lambda values.

- lambda.min.ratio:

  Ratio of min/max lambda.

- nfolds:

  Number of CV folds.

- cv.criteria:

  Cross-validation criteria.

- c_index_stratum:

  Stratum for C-index calculation.

- message:

  Logical. Print progress.

- seed:

  Seed for reproducibility.

- ...:

  Additional arguments for cv.coxkl_enet.

## Value

An object of class "bagging".

## Examples

``` r
if (FALSE) { # \dontrun{
data(ExampleData_highdim)
train_dat_highdim <- ExampleData_highdim$train
beta_external_highdim <- ExampleData_highdim$beta_external
etas <- generate_eta(method = "exponential", n = 10, max_eta = 100)

bagging.beta_fixed <- coxkl_enet_bagging(
  z = train_dat_highdim$z,
  delta = train_dat_highdim$status,
  time = train_dat_highdim$time,
  stratum = train_dat_highdim$stratum,
  beta = beta_external_highdim,
  etas = etas,
  B = 5,
  seed = 1
)
} # }
```
