# Example Data for Conditional Logistic Regression

A simulated dataset generated for 1:M matched case-control studies
(Conditional Logistic Regression, CLR). The data is organized into
matched sets (strata), with exactly one case (`y=1`) and \\m=4\\
controls (`y=0`) per set.

## Usage

``` r
data(ExampleData_cc)
```

## Format

A list containing the following elements:

- train:

  A list with components for training the CLR model:

  z

  :   Numeric matrix of covariates (dimension
      \\n\_{\mathrm{train}}\times 6\\) with columns named `Z1`–`Z6`.

  y

  :   Binary outcome vector (`1`=case, `0`=control).

  stratum

  :   Integer vector identifying the matched set for each observation
      (200 unique strata in `train`).

- test:

  A list with the same structure as `train`, used for external
  evaluation (500 unique strata in `test`).

- beta_external:

  Numeric vector (length 6) of CLR coefficients estimated on a separate
  external dataset using all `Z1`–`Z6`.

## Details

Data-generating mechanism:

- **Study design:** 1:4 matched case-control study (\\m=4\\ controls per
  case).

- **Covariates:** 6 variables (`Z1`–`Z6`) drawn from a correlated
  multivariate normal distribution.

- **True coefficients:** \\\beta = (1, -1, 1, -1, 1, -1)\\.

- **Set-specific effect:** A random stratum-specific intercept
  \\\theta_i \sim N(0, 0.5^2)\\ is added to the linear predictor; it is
  eliminated by CLR conditioning.

- **Outcome generation:** Within each stratum \\i\\, the single case
  (`y=1`) is selected with probability proportional to \\\exp(\theta_i +
  Z^\top \beta)\\.

- **External beta estimation:** `beta_external` is obtained by fitting
  `clogit` on a separate simulated dataset with a slightly different
  true coefficient vector \\\beta\_{\mathrm{ext}} = (0.8, -0.8, \dots)\\
  and correlation \\\rho = 0.3\\, using the `"breslow"` tie
  approximation.

## Examples

``` r
data(ExampleData_cc)

# Inspect the first few rows of the covariate matrix (z)
head(ExampleData_cc$train$z)
#>               Z1         Z2         Z3          Z4          Z5           Z6
#> [1,]  1.42087877  1.6308974  0.9538616  0.07376305  0.03428315  0.546671598
#> [2,]  1.82713082  0.8303178 -0.3442864 -0.48402535  0.48166198 -0.116491439
#> [3,]  1.22007941  0.7867172  1.6578626 -0.38287991 -2.76409690  1.818135337
#> [4,] -0.01727926 -0.5445138  0.2732903  0.83712396  0.64761711 -0.617661015
#> [5,]  1.02388718  0.8927722  0.8106354  1.17493045 -1.19407166 -0.009421722
#> [6,] -2.44226899  0.4824824 -0.1075286  0.34754847 -0.30607005  0.463093251

# Check the distribution of the binary outcome (y)
table(ExampleData_cc$train$y)
#> 
#>   0   1 
#> 800 200 
```
