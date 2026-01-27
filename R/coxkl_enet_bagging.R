#' Bagging for cv.coxkl_enet
#'
#' @param z Matrix of predictors (n x p).
#' @param delta Vector of event indicators.
#' @param time Vector of survival times.
#' @param stratum Vector indicating the stratum.
#' @param RS Vector of external risk scores. Resampled during bagging.
#' @param beta Vector of fixed external coefficients. Not resampled.
#' @param etas Sequence of eta values.
#' @param alpha Elastic net mixing parameter.
#' @param B Number of bootstrap replicates.
#' @param lambda Optional lambda sequence.
#' @param nlambda Number of lambda values.
#' @param lambda.min.ratio Ratio of min/max lambda.
#' @param nfolds Number of CV folds.
#' @param cv.criteria Cross-validation criteria.
#' @param c_index_stratum Stratum for C-index calculation.
#' @param message Logical. Print progress.
#' @param seed Seed for reproducibility.
#' @param ... Additional arguments for cv.coxkl_enet.
#'
#' @return An object of class "bagging".
#'
#' @examples
#' \dontrun{
#' data(ExampleData_highdim)
#' train_dat_highdim <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' etas <- generate_eta(method = "exponential", n = 10, max_eta = 100)
#'
#' bagging.beta_fixed <- coxkl_enet_bagging(
#'   z = train_dat_highdim$z,
#'   delta = train_dat_highdim$status,
#'   time = train_dat_highdim$time,
#'   stratum = train_dat_highdim$stratum,
#'   beta = beta_external_highdim,
#'   etas = etas,
#'   B = 5,
#'   seed = 1
#' )
#' }
#' @export
coxkl_enet_bagging <- function(z, delta, time, stratum = NULL, RS = NULL, beta = NULL,
                               etas, alpha = 1.0, B = 100, lambda = NULL, nlambda = 100,
                               lambda.min.ratio = ifelse(nrow(z) < ncol(z), 0.05, 1e-03),
                               nfolds = 5, cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
                               c_index_stratum = NULL,
                               message = FALSE, seed = NULL, ...) {

  cv.criteria <- match.arg(cv.criteria)

  z <- as.matrix(z)
  n <- nrow(z)
  p <- ncol(z)

  if (is.null(RS) && is.null(beta)) {
    stop("No external information is provided. Either RS or beta must be provided.")
  }

  if (!is.null(RS)) RS <- as.matrix(RS)

  if (is.null(stratum)) {
    stratum_full <- rep(1, n)
  } else {
    stratum_full <- stratum
  }

  if (!is.null(seed)) set.seed(seed)

  res_list <- vector("list", B)

  if (message) pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (i in seq_len(B)) {
    boot_idx <- sort(sample(seq_len(n), size = n, replace = TRUE))

    z_b       <- z[boot_idx, , drop = FALSE]
    delta_b   <- delta[boot_idx]
    time_b    <- time[boot_idx]
    stratum_b <- stratum_full[boot_idx]

    RS_b <- NULL
    beta_in <- NULL

    if (!is.null(RS)) {
      RS_b <- RS[boot_idx, , drop = FALSE]
    } else {
      beta_in <- beta
    }

    c_idx_strat_b <- NULL
    if (!is.null(c_index_stratum)) {
      c_idx_strat_b <- c_index_stratum[boot_idx]
    }

    fit_res <- tryCatch({
      cv.coxkl_enet(
        z = z_b,
        delta = delta_b,
        time = time_b,
        stratum = stratum_b,
        RS = RS_b,
        beta = beta_in,
        etas = etas,
        alpha = alpha,
        lambda = lambda,
        nlambda = nlambda,
        lambda.min.ratio = lambda.min.ratio,
        nfolds = nfolds,
        cv.criteria = cv.criteria,
        c_index_stratum = c_idx_strat_b,
        message = FALSE,
        ...
      )
    }, error = function(e) {
      return(NULL)
    })

    if (!is.null(fit_res)) {
      res_list[[i]] <- as.vector(fit_res$best$best_beta)
    } else {
      res_list[[i]] <- rep(NA, p)
    }

    if (message) setTxtProgressBar(pb, i)
  }

  if (message) close(pb)

  res_mat <- do.call(cbind, res_list)
  valid_cols <- !apply(res_mat, 2, function(x) any(is.na(x)))

  if (sum(valid_cols) < B) {
    warning(sprintf("Only %d out of %d bootstrap replicates converged.", sum(valid_cols), B))
    res_mat <- res_mat[, valid_cols, drop = FALSE]
  }

  bagged_beta <- rowMeans(res_mat)

  structure(
    list(
      best_beta = bagged_beta,
      all_betas = res_mat,
      B = B,
      seed = seed,
      valid_replicates = sum(valid_cols)
    ),
    class = "bagging"
  )
}
