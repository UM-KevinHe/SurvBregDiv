# Evaluate NCC Model: Loss and C-Index

Evaluate NCC Model: Loss and C-Index

## Usage

``` r
test_eval_ncc(z_ncc, case, set_id, betahat, criteria = c("loss", "CIndex"))
```

## Arguments

- z_ncc:

  Matrix of covariates for NCC data (rows = subjects).

- case:

  Integer or logical vector (0/1) indicating cases.

- set_id:

  Vector of matched set identifiers.

- betahat:

  Numeric vector of estimated coefficients.

- criteria:

  "loss" or "CIndex".

## Value

Numeric performance metric.
