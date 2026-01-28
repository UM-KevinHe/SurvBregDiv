#' Bagging for MDTL-Integrated Cox Elastic-Net Models
#'
#' Performs bootstrap aggregation (bagging) for the Mahalanobis-distance–based
#' transfer-learning Cox elastic-net model (\code{cv.cox_MDTL_enet}) by repeatedly
#' refitting the model on bootstrap resamples of the internal dataset and
#' averaging the resulting fitted coefficient vectors. This procedure reduces
#' sampling variability and improves robustness relative to a single data split.
#'
#' External information is supplied via a fixed coefficient vector (\code{beta})
#' and, optionally, a weighting matrix (\code{vcov}). Both represent external
#' prior information and are \strong{not} resampled across replicates.
#'
#' @param z Matrix of predictors of dimension \code{n x p}.
#' @param delta Event indicator vector.
#' @param time Survival time vector.
#' @param stratum Optional stratum indicator vector for stratified Cox modeling.
#' @param beta External coefficient vector of length \code{p}. Treated as fixed
#'   prior information and not resampled across bootstrap replicates.
#' @param vcov Optional weighting matrix (\code{p x p}) used in the Mahalanobis
#'   distance formulation.
#' @param etas Vector of \code{eta} values for transfer-learning shrinkage.
#' @param alpha Elastic-net mixing parameter between \code{0} and \code{1}.
#'   \code{alpha = 1} corresponds to lasso; \code{alpha = 0} to ridge. Default is \code{1.0}.
#' @param B Number of bootstrap replicates. Default is \code{100}.
#' @param lambda Optional user-specified \code{lambda} sequence.
#' @param nlambda Number of \code{lambda} values to generate if \code{lambda} is not supplied.
#' @param lambda.min.ratio Ratio of the smallest to the largest \code{lambda} when generating a sequence.
#' @param nfolds Number of folds for inner cross-validation via \code{cv.cox_MDTL_enet}.
#' @param cv.criteria Cross-validation criterion used for selecting the optimal
#'   \code{(eta, lambda)} pair.
#' @param c_index_stratum Optional stratum assignment for stratified C-index
#'   evaluation (may differ from model stratification).
#' @param message Logical indicating whether to print progress. Default is \code{FALSE}.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments passed to \code{cv.cox_MDTL_enet}.
#'
#' @return
#' An object of class \code{"cox_MDTL_bagging"} containing:
#' \itemize{
#'   \item \code{best_beta} — aggregated coefficient estimate obtained by averaging
#'     across valid bootstrap replicates.
#'   \item \code{all_betas} — matrix of dimension \code{p x B_valid} containing
#'     coefficient vectors from each successful bootstrap fit.
#'   \item \code{B} — total number of requested bootstrap replicates.
#'   \item \code{valid_replicates} — number of successful (non-error) fits contributing to aggregation.
#'   \item \code{seed} — seed used for reproducibility (if supplied).
#' }
#'
#' @examples
#' \dontrun{
#' data(ExampleData_highdim)
#' train_dat_highdim     <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#'
#' etas <- generate_eta(method = "exponential", n = 10, max_eta = 10)
#'
#' bag.out <- cox_MDTL_enet_bagging(
#'   z            = train_dat_highdim$z,
#'   delta        = train_dat_highdim$status,
#'   time         = train_dat_highdim$time,
#'   stratum      = train_dat_highdim$stratum,
#'   beta         = beta_external_highdim,
#'   vcov         = NULL,
#'   etas         = etas,
#'   alpha        = 0.5,
#'   B            = 5,
#'   cv.criteria  = "CIndex_pooled",
#'   message      = TRUE,
#'   seed         = 123
#' )
#' }
#'
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
