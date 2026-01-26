# Simulate 1:m matched case-control data

Internal function to simulate 1:m matched data. Each stratum contains
exactly one case. Probabilities are proportional to `exp(eta)`.

## Usage

``` r
sim.matched_cc(n_stratum, m, beta, rho = 0.8, stratum_sd = 0.5, seed = NULL)
```

## Arguments

- n_stratum:

  Number of matched strata (sets).

- m:

  Number of controls per case (\>=1).

- beta:

  Numeric vector of coefficients (length p).

- rho:

  Correlation parameter for Z within stratum.

- stratum_sd:

  SD of stratum random intercepts.

- seed:

  Optional RNG seed.

## Value

A list containing matched data frames and simulation metadata.
