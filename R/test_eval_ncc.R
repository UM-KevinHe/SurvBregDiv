#' Evaluate NCC Model: Loss, C-Index, and Brier Score
#'
#' @param z_ncc Matrix of covariates for NCC data (rows = subjects).
#' @param case Integer or logical vector (0/1) indicating cases.
#' @param set_id Vector of matched set identifiers.
#' @param betahat Numeric vector of estimated coefficients.
#' @param criteria "loss", "CIndex", or "Brier".
#'
#' @return Numeric performance metric.
#'
#' @keywords internal
#' @export
test_eval_ncc <- function(z_ncc,
                          case,
                          set_id,
                          betahat,
                          criteria = c("loss", "CIndex", "Brier")) {
  criteria <- match.arg(criteria)
  z_ncc <- as.matrix(z_ncc)
  betahat <- as.matrix(betahat)

  lp <- drop(z_ncc %*% betahat)
  case <- as.integer(case)
  set_id <- as.factor(set_id)

  levs <- levels(set_id)
  ord_list <- vector("list", length(levs))

  for (k in seq_along(levs)) {
    s <- levs[k]
    idx <- which(set_id == s)
    yk <- case[idx]

    idx_case <- idx[yk == 1L]
    idx_ctrl <- idx[yk == 0L]

    if (length(idx_case) != 1L) {
      stop("Each matched set must contain exactly one case.")
    }

    ord_list[[k]] <- c(idx_case, idx_ctrl)
  }

  ord <- unlist(ord_list, use.names = FALSE)

  lp_ord    <- lp[ord]
  case_ord  <- case[ord]
  set_ord   <- set_id[ord]

  n_each_stratum <- as.numeric(table(set_ord))

  if (criteria == "loss") {
    delta <- case_ord
    loglik <- pl_cal_theta(
      lp = lp_ord,
      delta = delta,
      n_each_stratum = n_each_stratum
    )
    n_sets <- length(n_each_stratum)
    return(as.numeric(-2 * loglik / n_sets))
  }

  if (criteria == "CIndex") {
    auc_obj <- auc_ncc(
      y = case_ord,
      score = lp_ord,
      set_id = set_ord
    )
    return(as.numeric(auc_obj$auc))
  }

  if (criteria == "Brier") {
    p_hat <- prob_in_set(lp_ord, set_ord)
    return(mean((case_ord - p_hat)^2))
  }

  stop("Unsupported criteria for NCC: ", criteria)
}


auc_one_set <- function(y, score) {
  y <- as.integer(y)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)

  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


auc_ncc <- function(y, score, set_id) {
  y <- as.integer(y)
  score <- as.numeric(score)
  set_id <- as.factor(set_id)

  levs <- levels(set_id)
  numer <- 0
  denom <- 0

  for (k in levs) {
    idx <- which(set_id == k)
    yk <- y[idx]
    sk <- score[idx]

    n1 <- sum(yk == 1L)
    n0 <- sum(yk == 0L)
    if (n1 == 0L || n0 == 0L) next

    auc_k <- auc_one_set(yk, sk)
    if (is.na(auc_k)) next

    N_pairs <- n1 * n0
    numer <- numer + auc_k * N_pairs
    denom <- denom + N_pairs
  }

  if (denom == 0) {
    return(list(numer = NA_real_, denom = 0, auc = NA_real_))
  }

  list(
    numer = numer,
    denom = denom,
    auc   = numer / denom
  )
}


prob_in_set <- function(lp, set_id) {
  lp <- as.numeric(lp)
  set_id <- as.factor(set_id)

  levs <- levels(set_id)
  p <- numeric(length(lp))

  for (s in levs) {
    idx <- which(set_id == s)
    lp_s <- lp[idx]
    a <- max(lp_s)
    w <- exp(lp_s - a)
    p[idx] <- w / sum(w)
  }

  p
}
