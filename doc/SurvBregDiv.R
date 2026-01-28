## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(fig.alt = "Plot generated in survkl vignette")
library(dplyr)

knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 4,
  dpi = 120
)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("SurvBregDiv")

## ----eval=FALSE---------------------------------------------------------------
# require(devtools)
# require(remotes)
# remotes::install_github("UM-KevinHe/SurvBregDiv", ref = "main")

## -----------------------------------------------------------------------------
library(SurvBregDiv)

## -----------------------------------------------------------------------------
data(ExampleData_lowdim)

train  <- ExampleData_lowdim$train
test   <- ExampleData_lowdim$test

z      <- train$z
delta  <- train$status
time   <- train$time
strat  <- train$stratum

## -----------------------------------------------------------------------------
beta_ext <- ExampleData_lowdim$beta_external_fair

## -----------------------------------------------------------------------------
eta_list <- generate_eta(method = "exponential", n = 50, max_eta = 10)

## -----------------------------------------------------------------------------
coxkl_est <- coxkl(
  z = z,
  delta = delta,
  time = time,
  stratum = strat,
  beta = beta_ext,
  etas = eta_list
)

## -----------------------------------------------------------------------------
cox_MDTL_est <- cox_MDTL(
  z = z,
  delta = delta,
  time = time,
  stratum = strat,
  beta = beta_ext,
  vcov = NULL,
  etas = eta_list
)

## -----------------------------------------------------------------------------
RS_ext <- as.matrix(z) %*% as.matrix(beta_ext)

coxkl_est.RS <- coxkl(
  z = z,
  delta = delta,
  time = time,
  stratum = strat,
  RS = RS_ext,
  etas = eta_list
)

## -----------------------------------------------------------------------------
time_ties <- round(time, 2)   # Rounding time introduces ties for demonstration

coxkl_ties_est <- coxkl_ties(
  z = z,
  delta = delta,
  time = time_ties,
  stratum = strat,
  beta = beta_ext,
  etas = eta_list,
  ties = "breslow"
)

## -----------------------------------------------------------------------------
cv.coxkl_est <- cv.coxkl(
  z = z,
  delta = delta,
  time = time,
  stratum = strat,
  beta = beta_ext,
  etas = eta_list,
  nfolds = 5,
  criteria = "V&VH",
  seed = 1)

## -----------------------------------------------------------------------------
plot(
  cox_MDTL_est,
  test_z       = test$z,
  test_time    = test$time,
  test_delta   = test$status,
  test_stratum = test$stratum,
  criteria     = "loss"
) 

## -----------------------------------------------------------------------------
cv.plot(cv.coxkl_est)

## -----------------------------------------------------------------------------
data(ExampleData_highdim)

train_hd  <- ExampleData_highdim$train
test_hd   <- ExampleData_highdim$test

z_hd      <- train_hd$z
delta_hd  <- train_hd$status
time_hd   <- train_hd$time
strat_hd  <- train_hd$stratum

beta_external_hd <- ExampleData_highdim$beta_external

## -----------------------------------------------------------------------------
model_ridge <- coxkl_ridge(
  z        = z_hd,
  delta    = delta_hd,
  time     = time_hd,
  stratum  = strat_hd,
  beta     = beta_external_hd,
  eta      = 1                 
)

## -----------------------------------------------------------------------------
model_enet <- coxkl_enet(
  z        = z_hd,
  delta    = delta_hd,
  time     = time_hd,
  stratum  = strat_hd,
  beta     = beta_external_hd,   
  eta      = 1,                 
  alpha    = 1                  # LASSO penalty
)

## -----------------------------------------------------------------------------
eta_grid_hd <- generate_eta(method = "exponential",
                            n = 50,
                            max_eta = 100)

cv_enet_hd <- cv.coxkl_enet(
  z = z_hd,
  delta = delta_hd,
  time = time_hd,
  stratum = strat_hd,
  beta = beta_external_hd,
  etas = eta_grid_hd,
  alpha = 1, # LASSO
  nfolds = 5,
  cv.criteria = "V&VH",
  seed = 1
)

## -----------------------------------------------------------------------------
cv_enet_hd$integrated_stat.best_per_eta

## -----------------------------------------------------------------------------
plot(
  model_enet,
  test_z       = test_hd$z,
  test_time    = test_hd$time,
  test_delta   = test_hd$status,
  test_stratum = test_hd$stratum,
  criteria     = "loss"
)

## -----------------------------------------------------------------------------
cv.plot(cv_enet_hd)

## -----------------------------------------------------------------------------
coxkl.StabSelect <- coxkl_enet.StabSelect(
  z = z_hd,
  delta = delta_hd,
  time = time_hd,
  stratum = strat_hd,
  beta = beta_external_hd,
  etas = eta_list,
  cv.criteria = "CIndex_pooled",
  B = 20
)

## -----------------------------------------------------------------------------
plot(coxkl.StabSelect, threshold = 0.6) 

## ----eval=FALSE, message=FALSE------------------------------------------------
# bagging.coxkl <- coxkl_enet_bagging(
#   z = z_hd,
#   delta = delta_hd,
#   time = time_hd,
#   stratum = strat_hd,
#   beta = beta_external_hd,
#   etas = eta_list,
#   B = 5,
#   seed = 1
# )

## -----------------------------------------------------------------------------
data(ExampleData_cc)

train.cc  <- ExampleData_cc$train
test.cc   <- ExampleData_cc$test

z.cc      <- train.cc$z
y.cc      <- train.cc$y
set.cc    <- train.cc$stratum

beta_ext.cc <- ExampleData_cc$beta_external

## -----------------------------------------------------------------------------
eta_list <- generate_eta(method = "exponential", n = 50, max_eta = 5)

clogitkl.fit_breslow <- clogitkl(y = y.cc, z = z.cc, stratum = set.cc, 
                                 eta = eta_list, beta = beta_ext.cc,
                                 method = "breslow")

## -----------------------------------------------------------------------------
plot(
  clogitkl.fit_breslow,
  test_z       = test.cc$z,
  test_delta   = test.cc$y,
  test_stratum = test.cc$stratum,
  criteria     = "loss"
)

## -----------------------------------------------------------------------------
cv.clogitkl.fit_breslow <- cv.clogitkl(
  y        = y.cc,
  z        = z.cc,
  stratum  = set.cc,
  beta     = beta_ext.cc,
  etas     = eta_list,
  method   = "exact",
  nfolds   = 5,
  criteria = "loss"
)

## -----------------------------------------------------------------------------
cv.plot(cv.clogitkl.fit_breslow)

## -----------------------------------------------------------------------------
data(ExampleData_cc_highdim)

train.cc_hd  <- ExampleData_cc_highdim$train
test.cc_hd   <- ExampleData_cc_highdim$test

z.cc_hd      <- train.cc_hd$z
y.cc_hd      <- train.cc_hd$y
set.cc_hd    <- train.cc_hd$stratum

beta_ext.cc_hd <- ExampleData_cc_highdim$beta_external

## -----------------------------------------------------------------------------
clogitkl_enet_fit <- clogitkl_enet(
  y       = y.cc_hd,
  z       = z.cc_hd,
  stratum = set.cc_hd,
  beta    = beta_ext.cc_hd,
  eta     = 0
)

## -----------------------------------------------------------------------------
plot(
  clogitkl_enet_fit,
  test_z       = test.cc_hd$z,
  test_delta   = test.cc_hd$y,
  test_stratum = test.cc_hd$stratum,
  criteria     = "loss"
)

## -----------------------------------------------------------------------------
eta_list <- generate_eta(method = "exponential", n = 50, max_eta = 5)

cv.clogitkl_enet_fit <- cv.clogitkl_enet(
  y        = y.cc_hd,
  z        = z.cc_hd,
  stratum  = set.cc_hd,
  beta     = beta_ext.cc_hd,
  etas     = eta_list,
  alpha    = 1,
  nfolds   = 5,
  criteria = "loss"
)


## -----------------------------------------------------------------------------
cv.plot(cv.clogitkl_enet_fit)

