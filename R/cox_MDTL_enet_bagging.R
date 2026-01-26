#' Bagging for cv.cox_MDTL_enet
#'
#' @description
#' Implements bootstrap aggregation (bagging) for the \code{cv.cox_MDTL_enet} model.
#' It generates \code{B} bootstrap samples from the internal data, fits the Cross-Validated
#' Cox MDTL Elastic Net model on each sample, and aggregates the resulting coefficients.
#'
#' @param z Matrix of predictors (n x p).
#' @param delta Vector of event indicators.
#' @param time Vector of survival times.
#' @param stratum Vector indicating the stratum.
#' @param beta Vector of fixed external coefficients (length p). This is prior information and is \strong{not} resampled.
#' @param vcov Optional weighting matrix (p x p).
#' @param etas Sequence of eta values (transfer learning weights) to tune.
#' @param alpha The Elastic Net mixing parameter, with \eqn{0 \le \alpha \le 1}.
#'   \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty. Default is 1.0.
#' @param B Integer. Number of bootstrap replicates. Default is 100.
#' @param lambda Optional user-supplied lambda sequence.
#' @param nlambda Number of lambda values.
#' @param lambda.min.ratio Ratio of min/max lambda.
#' @param nfolds Number of CV folds for the inner cross-validation.
#' @param cv.criteria Cross-validation criteria.
#' @param c_index_stratum Stratum vector for C-index calculation (if different from model stratum).
#' @param message Logical. If TRUE, shows a progress bar.
#' @param seed Integer. Seed for reproducibility.
#' @param ... Additional arguments passed to \code{cv.cox_MDTL_enet}.
#'
#' @return An object of class \code{"cox_MDTL_bagging"} containing:
#' \itemize{
#'   \item \code{best_beta}: The averaged coefficient vector across all valid bootstrap replicates.
#'   \item \code{all_betas}: A matrix (p x B) of coefficients from each bootstrap replicate.
#'   \item \code{B}: Number of requested replicates.
#'   \item \code{valid_replicates}: Number of replicates that successfully converged.
#'   \item \code{seed}: The seed used.
#' }
#'
#' @examples
#' \donttest{
#' data(ExampleData_highdim)
#' train_dat_highdim <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' 
#' etas <- generate_eta(method = "exponential", n = 10, max_eta = 10)
#' bagging_res <- cox_MDTL_enet_bagging(
#'   z = train_dat_highdim$z,
#'   delta = train_dat_highdim$status,
#'   time = train_dat_highdim$time,
#'   stratum = train_dat_highdim$stratum,
#'   beta = beta_external_highdim,
#'   vcov = NULL,
#'   etas = etas,
#'   alpha = 0.5, # Elastic Net mixing
#'   B = 5,
#'   cv.criteria = "CIndex_pooled",
#'   message = TRUE,
#'   seed = 123
#' )
#' }
#' @export
cox_MDTL_enet_bagging <- function(z, delta, time, stratum = NULL, beta = NULL, vcov = NULL,
                                  etas, alpha = 1.0, B = 100, lambda = NULL, nlambda = 100,
                                  lambda.min.ratio = ifelse(nrow(z) < ncol(z), 0.01, 1e-04),
                                  nfolds = 5,
                                  cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
                                  c_index_stratum = NULL,
                                  message = FALSE, seed = NULL, ...) {
  
  cv.criteria <- match.arg(cv.criteria)
  
  z <- as.matrix(z)
  n <- nrow(z)
  p <- ncol(z)
  
  # Input checks specific to MDTL
  if (is.null(beta)) {
    stop("External beta must be provided for Cox MDTL.")
  }
  if (length(beta) != p) {
    stop("Length of external beta must match number of columns in z.")
  }
  
  if (is.null(stratum)) {
    stratum_full <- rep(1, n)
  } else {
    stratum_full <- stratum
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  res_list <- vector("list", B)
  
  if (message) {
    cat("Starting Bagging (B =", B, ") for cv.cox_MDTL_enet:\n")
    pb <- txtProgressBar(min = 0, max = B, style = 3)
  }
  
  for (i in seq_len(B)) {
    # 1. Bootstrap sampling of the INTERNAL data
    boot_idx <- sort(sample(seq_len(n), size = n, replace = TRUE))
    
    z_b       <- z[boot_idx, , drop = FALSE]
    delta_b   <- delta[boot_idx]
    time_b    <- time[boot_idx]
    stratum_b <- stratum_full[boot_idx]
    
    # Note: 'beta' and 'vcov' are EXTERNAL information, so they are fixed and NOT resampled.
    
    c_idx_strat_b <- NULL
    if (!is.null(c_index_stratum)) {
      c_idx_strat_b <- c_index_stratum[boot_idx]
    }
    
    # 2. Fit cv.cox_MDTL_enet on bootstrapped data
    fit_res <- tryCatch({
      cv.cox_MDTL_enet(
        z = z_b,
        delta = delta_b,
        time = time_b,
        stratum = stratum_b,
        beta = beta,      # Fixed external beta
        vcov = vcov,      # Fixed external vcov
        etas = etas,
        alpha = alpha,    # Passed alpha
        lambda = lambda,
        nlambda = nlambda,
        lambda.min.ratio = lambda.min.ratio,
        nfolds = nfolds,
        cv.criteria = cv.criteria,
        c_index_stratum = c_idx_strat_b,
        message = FALSE,  # Suppress inner loop messages
        ...
      )
    }, error = function(e) {
      # Handle convergence failures gracefully
      return(NULL)
    })
    
    # 3. Store result
    if (!is.null(fit_res)) {
      # Extract the beta corresponding to the best (eta, lambda) combination
      res_list[[i]] <- as.vector(fit_res$best$best_beta)
    } else {
      res_list[[i]] <- rep(NA, p)
    }
    
    if (message) setTxtProgressBar(pb, i)
  }
  
  if (message) close(pb)
  
  # 4. Aggregate results
  res_mat <- do.call(cbind, res_list)
  
  # Check for failed runs (NA columns)
  valid_cols <- !apply(res_mat, 2, function(x) any(is.na(x)))
  n_valid <- sum(valid_cols)
  
  if (n_valid < B) {
    warning(sprintf("Only %d out of %d bootstrap replicates converged.", n_valid, B))
    res_mat <- res_mat[, valid_cols, drop = FALSE]
  }
  
  if (n_valid == 0) {
    stop("All bootstrap replicates failed.")
  }
  
  bagged_beta <- rowMeans(res_mat)
  
  structure(
    list(
      best_beta = bagged_beta,
      all_betas = res_mat,
      B = B,
      seed = seed,
      valid_replicates = n_valid
    ),
    class = "cox_MDTL_bagging"
  )
}