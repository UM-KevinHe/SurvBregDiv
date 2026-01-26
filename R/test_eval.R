#' Evaluate Survival Model Performance
#'
#' Computes predictive performance metrics for stratified or unstratified Cox models.
#' Supports Loss, C-index, Integrated Brier Score (IBS), and Time-Dependent AUC (tdAUC).
#'
#' @param test_z Matrix of predictors for the test set.
#' @param test_delta Numeric vector of event indicators (1 for event, 0 for censored).
#' @param test_time Numeric vector of observed times.
#' @param betahat Numeric vector of estimated coefficients.
#' @param test_stratum Vector indicating strata for test subjects. Defaults to NULL (single stratum).
#' @param train_baseline_obj A list containing the baseline hazard function (typically from \code{get_baseline_hazard}).
#' Required only when \code{criteria = "IBS"}.
#' @param criteria Metric to calculate: "loss" (Log-Partial Likelihood), "CIndex" (Concordance Index),
#' "IBS" (Integrated Brier Score), or "tdAUC" (Integrated Time-Dependent AUC).
#'
#' @details
#' For "IBS", the function predicts survival probabilities and converts them to risk (1 - S).
#' If \code{riskRegression} fails to provide a pre-computed IBS, the function manually integrates
#' the Brier score using the trapezoidal rule.
#'
#' @return A numeric value representing the performance metric.
#' Returns \code{NA} if the metric cannot be computed (e.g., no events in test set).
#'
#' @importFrom riskRegression Score
#' @importFrom survival Surv
#' @importFrom stats approxfun plnorm plogis pnorm qnorm rbinom relevel rexp rlnorm rnorm rpois runif rweibull
#' @importFrom utils head tail
#'
#' @keywords internal
#' @export
test_eval <- function(test_z, test_delta, test_time,
                      betahat, test_stratum = NULL,
                      train_baseline_obj = NULL,
                      criteria = c("loss", "CIndex", "IBS", "tdAUC")) {

  criteria <- match.arg(criteria)
  test_RS <- as.vector(as.matrix(test_z) %*% as.matrix(betahat))

  d_test <- data.frame(time = as.numeric(test_time), status = as.numeric(test_delta))
  n <- nrow(d_test)
  if (is.null(test_stratum)) test_stratum <- rep(1, n)

  if (criteria == "loss") {
    ord <- order(test_stratum, d_test$time)
    return(-2 * pl_cal_theta(test_RS[ord], d_test$status[ord], as.numeric(table(test_stratum))) / n)
  }

  if (criteria == "CIndex") {
    return(c_stat_stratcox(d_test$time, test_RS, test_stratum, d_test$status)$c_statistic)
  }

  eval_times <- sort(unique(d_test$time[d_test$status == 1]))

  # if (criteria == "IBS") {
  #   S_mat <- predict_surv_prob(
  #     test_RS = test_RS,
  #     eval_times = eval_times,
  #     train_baseline_obj = train_baseline_obj,
  #     test_stratum = test_stratum
  #   )
  #
  #   risk_mat <- 1 - S_mat
  #
  #   score_res <- riskRegression::Score(
  #     object = list("MyModel" = risk_mat),
  #     formula = Surv(time, status) ~ 1,
  #     data = d_test,
  #     times = eval_times,
  #     metrics = "brier",
  #     summary = "ibs"
  #   )
  #
  #   brier_summary <- as.data.frame(score_res$Brier$summary)
  #
  #   if (nrow(brier_summary) > 0 && "IBS" %in% colnames(brier_summary)) {
  #     ibs_val <- brier_summary[brier_summary$model == "MyModel", "IBS"]
  #   } else {
  #     brier_scores <- as.data.frame(score_res$Brier$score)
  #     brier_scores <- brier_scores[brier_scores$model == "MyModel", ]
  #     x <- brier_scores$times
  #     y <- brier_scores$Brier
  #     ibs_val <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
  #   }
  #
  #   return(as.numeric(ibs_val))
  # }


  if (criteria == "IBS") {
    S_mat <- predict_surv_prob(
      test_RS = test_RS,
      eval_times = eval_times,
      train_baseline_obj = train_baseline_obj,
      test_stratum = test_stratum
    )
    risk_mat <- 1 - S_mat

    unique_strata <- sort(unique(test_stratum))
    n_total <- nrow(d_test)
    ibs_vals <- numeric(length(unique_strata))
    weights <- numeric(length(unique_strata))

    for (k in seq_along(unique_strata)) {
      s <- unique_strata[k]
      idx <- which(test_stratum == s)

      if (length(idx) < 2L || all(d_test$status[idx] == 0)) {
        ibs_vals[k] <- NA_real_
        weights[k]  <- 0
        next
      }

      d_s <- d_test[idx, , drop = FALSE]
      risk_s <- risk_mat[idx, , drop = FALSE]

      score_res_s <- riskRegression::Score(
        object  = list("MyModel" = risk_s),
        formula = survival::Surv(time, status) ~ 1,
        data    = d_s,
        times   = eval_times,
        metrics = "brier",
        summary = "ibs"
      )

      brier_summary_s <- as.data.frame(score_res_s$Brier$summary)

      if (nrow(brier_summary_s) > 0 && "IBS" %in% colnames(brier_summary_s)) {
        ibs_vals[k] <- brier_summary_s[
          brier_summary_s$model == "MyModel",
          "IBS"
        ]
        weights[k] <- length(idx)
      } else {
        ibs_vals[k] <- NA_real_
        weights[k]  <- 0
      }
    }

    if (sum(weights) == 0) return(NA_real_)

    ibs_val <- sum(ibs_vals * weights, na.rm = TRUE) / sum(weights)
    return(as.numeric(ibs_val))
  }

  # if (criteria == "tdAUC") {
  #   score_res <- riskRegression::Score(
  #     object = list("MyModel" = test_RS),
  #     formula = Surv(time, status) ~ 1,
  #     data = d_test,
  #     times = eval_times,
  #     metrics = "auc"
  #   )
  #
  #   auc_df <- as.data.frame(score_res$AUC$score)
  #   auc_df <- auc_df[auc_df$model == "MyModel", ]
  #
  #   x <- auc_df$time
  #   y <- auc_df$AUC
  #   valid <- !is.na(y)
  #   if(sum(valid) < 2) return(mean(y, na.rm = TRUE))
  #
  #   iauc <- sum(diff(x[valid]) * (head(y[valid], -1) + tail(y[valid], -1)) / 2) /
  #     (max(x[valid]) - min(x[valid]))
  #   return(iauc)
  # }


  if (criteria == "tdAUC") {
    # Global event times
    eval_times <- sort(unique(d_test$time[d_test$status == 1]))
    if (length(eval_times) < 2L) {
      return(NA_real_)
    }

    unique_strata <- sort(unique(test_stratum))
    n_times <- length(eval_times)
    n_strata <- length(unique_strata)

    # dN_s(t): rows = strata, cols = times
    dN_mat <- matrix(0, nrow = n_strata, ncol = n_times)
    rownames(dN_mat) <- as.character(unique_strata)

    for (k in seq_along(unique_strata)) {
      s <- unique_strata[k]
      idx <- which(test_stratum == s)
      time_s <- d_test$time[idx]
      status_s <- d_test$status[idx]

      for (j in seq_along(eval_times)) {
        t <- eval_times[j]
        dN_mat[k, j] <- sum(time_s == t & status_s == 1)
      }
    }

    # AUC_s(t): same dimension as dN_mat
    AUC_mat <- matrix(NA_real_, nrow = n_strata, ncol = n_times)
    rownames(AUC_mat) <- as.character(unique_strata)

    for (k in seq_along(unique_strata)) {
      s <- unique_strata[k]
      idx <- which(test_stratum == s)

      # Skip strata with fewer than 2 subjects or no events
      if (length(idx) < 2L || all(d_test$status[idx] == 0)) {
        next
      }

      d_s <- d_test[idx, , drop = FALSE]
      rs  <- test_RS[idx]

      score_res_s <- riskRegression::Score(
        object  = list("MyModel" = rs),
        formula = survival::Surv(time, status) ~ 1,
        data    = d_s,
        times   = eval_times,
        metrics = "auc"
      )

      auc_df_s <- as.data.frame(score_res_s$AUC$score)
      auc_df_s <- auc_df_s[auc_df_s$model == "MyModel", ]

      # Align AUC_s(t) with eval_times
      AUC_vec <- rep(NA_real_, n_times)
      match_idx <- match(eval_times, auc_df_s$time)
      ok <- !is.na(match_idx)
      AUC_vec[ok] <- auc_df_s$AUC[match_idx[ok]]

      AUC_mat[k, ] <- AUC_vec
    }

    # Aggregate across strata using dN_s(t) weights
    numer_t <- colSums(dN_mat * AUC_mat, na.rm = TRUE)
    denom_t <- colSums(dN_mat, na.rm = TRUE)

    AUC_t <- ifelse(denom_t > 0, numer_t / denom_t, NA_real_)

    valid <- !is.na(AUC_t)
    if (sum(valid) < 2L) {
      return(mean(AUC_t, na.rm = TRUE))
    }

    x <- eval_times[valid]
    y <- AUC_t[valid]

    iauc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
    return(iauc)
  }

}


predict_surv_prob <- function(test_RS, eval_times, train_baseline_obj, test_stratum = NULL) {
  n <- length(test_RS)
  if (is.null(test_stratum)) test_stratum <- rep(1, n)

  S_pred <- matrix(0, nrow = n, ncol = length(eval_times))
  unique_strata <- unique(test_stratum)

  for (s in unique_strata) {
    idx <- which(test_stratum == s)
    Lambda0_t <- train_baseline_obj$predict_baseline(eval_times, strat_id = s)
    S_pred[idx, ] <- exp(-exp(test_RS[idx]) %o% Lambda0_t)
  }
  return(S_pred)
}

get_baseline_hazard <- function(z, delta, time, beta, stratum = NULL) {
  lp <- drop(as.matrix(z) %*% as.matrix(beta))

  if (is.null(stratum)) {
    df <- data.frame(time = time, status = delta)
    fit_fixed <- survival::coxph(Surv(time, status) ~ offset(lp), data = df)
    bh <- survival::basehaz(fit_fixed, centered = FALSE)

    fun <- approxfun(
      x = bh$time,
      y = bh$hazard,
      method = "constant",
      yleft = 0,
      rule = 2
    )

    predict_baseline <- function(times, strat_id = NULL) {
      fun(times)
    }

    return(list(predict_baseline = predict_baseline))
  }

  str_raw <- stratum
  u <- sort(unique(str_raw))

  baseline_funs <- vector("list", length(u))
  names(baseline_funs) <- as.character(u)

  for (s in u) {
    idx <- which(str_raw == s)
    df_s <- data.frame(time = time[idx], status = delta[idx])

    fit_s <- survival::coxph(
      Surv(time, status) ~ offset(lp[idx]),
      data = df_s
    )
    bh_s <- survival::basehaz(fit_s, centered = FALSE)

    baseline_funs[[as.character(s)]] <- approxfun(
      x = bh_s$time,
      y = bh_s$hazard,
      method = "constant",
      yleft = 0,
      rule = 2
    )
  }

  predict_baseline <- function(times, strat_id) {
    sid <- as.character(strat_id)
    if (!sid %in% names(baseline_funs)) {
      stop("Stratum ", strat_id, " not present in training data.")
    }
    baseline_funs[[sid]](times)
  }

  list(predict_baseline = predict_baseline)
}































