#' Multi-Source Integration for KL-Integrated Cox Elastic-Net Models
#'
#' Fits multiple KL-integrated Cox elastic-net models on the full data using multiple
#' external sources, and combines the fitted coefficient vectors across sources to
#' produce a single aggregated estimate.
#'
#' Unlike \code{coxkl_enet_bagging()}, this function does not bootstrap the data.
#' Instead, it runs \code{cv.coxkl_enet()} once per external source on the full dataset.
#' The resulting coefficient vectors are then aggregated (by default, averaged) to
#' obtain a combined estimate.
#'
#' @param z Matrix/data.frame of predictors of dimension \code{n x p}.
#' @param delta Event indicator vector.
#' @param time Survival time vector.
#' @param stratum Optional stratum indicator vector for stratified Cox models.
#' @param beta_list A list of external coefficient vectors. Each element must have
#'   length \code{p}. If provided, \code{RS_list} should be \code{NULL}.
#' @param RS_list Optional list of external risk score vectors/matrices. Each element
#'   should be conformable with \code{n}. If provided, \code{beta_list} is ignored.
#' @param etas Vector of \code{eta} values for transfer-learning shrinkage.
#' @param combine How to combine coefficients across sources. Either \code{"mean"} (default)
#'   or \code{"median"}.
#' @param message Logical indicating whether to print progress.
#' @param seed Optional seed for reproducibility (passed to each CV run with an offset).
#' @param ... Additional arguments passed to \code{cv.coxkl_enet()} (e.g., \code{alpha},
#'   \code{lambda}, \code{nlambda}, \code{lambda.min.ratio}, \code{nfolds},
#'   \code{cv.criteria}, \code{c_index_stratum}, etc.).
#'
#' @return
#' An object of class \code{"coxkl_enet.multi"}, which is a list containing:
#' \itemize{
#'   \item \code{best_beta} — combined coefficient estimate across sources.
#'   \item \code{all_betas} — matrix of dimension \code{p x K_valid} of coefficient vectors
#'     from each successful fit.
#'   \item \code{K} — total number of external sources provided.
#'   \item \code{valid_sources} — number of successful (non-error) fits used in aggregation.
#'   \item \code{combine} — combination rule used.
#'   \item \code{seed} — seed used (if any).
#' }
#'
#' @export
coxkl_enet.multi <- function(
    z, delta, time, stratum = NULL,
    beta_list = NULL, RS_list = NULL,
    etas,
    combine = c("mean", "median"),
    message = FALSE,
    seed = NULL,
    ...
) {
  combine <- match.arg(combine)

  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z)) stop("z must be a matrix or a data.frame convertible to a matrix.", call. = FALSE)

  storage.mode(z) <- "double"
  if (!is.numeric(z)) stop("z must be numeric after conversion.", call. = FALSE)
  if (anyNA(z)) stop("z contains NA values. Please impute or remove missing values before fitting.", call. = FALSE)

  n <- nrow(z)
  p <- ncol(z)

  if (length(delta) != n) stop("delta length must match nrow(z).", call. = FALSE)
  if (length(time) != n) stop("time length must match nrow(z).", call. = FALSE)
  if (!is.null(stratum) && length(stratum) != n) stop("stratum length must match nrow(z).", call. = FALSE)

  if (!is.null(RS_list)) {
    if (!is.list(RS_list) || length(RS_list) < 1) stop("RS_list must be a non-empty list.", call. = FALSE)
    K <- length(RS_list)
  } else {
    if (!is.list(beta_list) || length(beta_list) < 1) stop("beta_list must be a non-empty list.", call. = FALSE)
    K <- length(beta_list)
    for (k in seq_len(K)) {
      if (length(beta_list[[k]]) != p) stop("beta_list element dimension mismatch with z.", call. = FALSE)
    }
  }

  stratum_full <- if (is.null(stratum)) rep(1, n) else stratum

  if (message) {
    cat("Fitting integrated models over external sources:\n")
    pb <- txtProgressBar(min = 0, max = K, style = 3, width = 30)
  }

  res_list <- vector("list", K)

  for (k in seq_len(K)) {
    RS_k <- NULL
    beta_k <- NULL

    if (!is.null(RS_list)) {
      RS_k <- as.matrix(RS_list[[k]])
      if (nrow(RS_k) != n) stop("Each RS_list element must have nrow equal to nrow(z).", call. = FALSE)
    } else {
      beta_k <- beta_list[[k]]
    }

    fit_res <- tryCatch({
      cv.coxkl_enet(
        z = z,
        delta = delta,
        time = time,
        stratum = stratum_full,
        RS = RS_k,
        beta = beta_k,
        etas = etas,
        message = FALSE,
        seed = if (is.null(seed)) NULL else (seed + k),
        ...
      )
    }, error = function(e) {
      return(NULL)
    })

    if (!is.null(fit_res)) {
      res_list[[k]] <- as.vector(fit_res$best$best_beta)
    } else {
      res_list[[k]] <- rep(NA_real_, p)
    }

    if (message) setTxtProgressBar(pb, k)
  }

  if (message) close(pb)

  res_mat <- do.call(cbind, res_list)
  valid_cols <- !apply(res_mat, 2, function(x) any(is.na(x)))

  if (sum(valid_cols) < K) {
    warning(sprintf("Only %d out of %d external sources produced successful fits.", sum(valid_cols), K))
    res_mat <- res_mat[, valid_cols, drop = FALSE]
  }

  if (ncol(res_mat) == 0) stop("No successful fits were obtained.", call. = FALSE)

  if (combine == "mean") {
    combined_beta <- rowMeans(res_mat)
  } else {
    combined_beta <- apply(res_mat, 1, stats::median)
  }

  structure(
    list(
      best_beta = combined_beta,
      all_betas = res_mat,
      K = K,
      seed = seed,
      valid_sources = sum(valid_cols),
      combine = combine
    ),
    class = "coxkl_enet.multi"
  )
}
