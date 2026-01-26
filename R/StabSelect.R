#' Stability Selection for Cox-KL Enet Model
#'
#' @inheritParams cv.coxkl_enet
#' @param B Integer. Number of bootstrap/subsampling replicates. Default is 50.
#' @param fraction_sample Numeric. Fraction of data to use for subsampling in each replicate. Default is 0.5.
#'
#' @return An object of class \code{StabSelect}, which is a list containing:
#' \item{stability_path}{A matrix (n_vars x n_lambda) of selection probabilities.}
#' \item{lambda}{The global lambda sequence used.}
#'
#' @examples
#' \dontrun{
#' data(ExampleData_highdim)
# train_dat_highdim <- ExampleData_highdim$train
# beta_external_highdim <- ExampleData_highdim$beta_external
# 
# eta_list <- generate_eta(method = "exponential", n = 10, max_eta = 50)
# 
# coxkl.StabSelect <- coxkl_enet.StabSelect(
#   z = train_dat_highdim$z,
#   delta = train_dat_highdim$status,
#   time = train_dat_highdim$time,
#   stratum = train_dat_highdim$stratum,
#   beta = beta_external_highdim,
#   etas = eta_list,
#   alpha = 1,
#   cv.criteria = "CIndex_pooled",
#   B = 20,
#   message = TRUE
# )
#'
#' # Plot with different thresholds without re-running
#' plot(coxkl.StabSelect, threshold = 0.6)
#' plot(coxkl.StabSelect, threshold = 0.8)
#' }
#'
#' @export
coxkl_enet.StabSelect <- function(z, delta, time, stratum = NULL, RS = NULL, beta = NULL,
                                  etas, alpha = 1.0,
                                  lambda = NULL, nlambda = 100, lambda.min.ratio = 0.1,
                                  nfolds = 5,
                                  cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
                                  c_index_stratum = NULL,
                                  message = FALSE, seed = NULL,
                                  B = 50,
                                  fraction_sample = 0.5,
                                  ...) {
  
  if (!is.null(seed)) set.seed(seed)
  cv.criteria <- match.arg(cv.criteria)
  
  z <- as.matrix(z)
  n_full <- nrow(z)
  p_vars <- ncol(z)
  var_names <- colnames(z)
  if(is.null(var_names)) var_names <- paste0("V", 1:p_vars)
  
  if (is.null(lambda.min.ratio)) {
    lambda.min.ratio <- ifelse(n_full < p_vars, 0.05, 1e-03)
  }
  
  if (is.null(lambda)) {
    if(message) message("Generating global lambda sequence using full data...")
    lambda_max_candidates <- numeric(length(etas))
    for(i in seq_along(etas)) {
      fit_tmp <- coxkl_enet(z = z, delta = delta, time = time, stratum = stratum,
                            RS = RS, beta = beta, eta = etas[i], alpha = alpha,
                            nlambda = 5, lambda.min.ratio = lambda.min.ratio, ...)
      if(!is.null(fit_tmp$lambda)) lambda_max_candidates[i] <- max(fit_tmp$lambda)
    }
    global_lambda_max <- 1.1 * max(lambda_max_candidates)
    lambda <- exp(seq(log(global_lambda_max), log(global_lambda_max * lambda.min.ratio), length.out = nlambda))
  }
  
  n_lambda_seq <- length(lambda)
  stability_counts <- matrix(0, nrow = p_vars, ncol = n_lambda_seq)
  rownames(stability_counts) <- var_names
  
  if(message) {
    message(sprintf("Starting Stability Selection with B=%d replicates...", B))
    pb <- txtProgressBar(min = 0, max = B, style = 3)
  }
  
  for (b in 1:B) {
    if(message) setTxtProgressBar(pb, b)
    
    sub_idx <- sort(sample(seq_len(n_full), size = floor(fraction_sample * n_full)))
    z_sub <- z[sub_idx, , drop = FALSE]
    delta_sub <- delta[sub_idx]
    time_sub <- time[sub_idx]
    stratum_sub <- if(!is.null(stratum)) stratum[sub_idx] else NULL
    RS_sub <- if(!is.null(RS)) RS[sub_idx] else NULL
    
    tryCatch({
      cv_fit <- cv.coxkl_enet(
        z = z_sub, delta = delta_sub, time = time_sub, stratum = stratum_sub,
        RS = RS_sub, beta = beta,
        etas = etas, alpha = alpha,
        lambda = lambda, 
        nfolds = nfolds,
        cv.criteria = cv.criteria,
        message = FALSE, 
        ...
      )
      
      best_eta <- cv_fit$best$best_eta
      
      fit_best <- coxkl_enet(
        z = z_sub, delta = delta_sub, time = time_sub, stratum = stratum_sub,
        RS = RS_sub, beta = beta,
        eta = best_eta,
        lambda = lambda,
        alpha = alpha,
        ...
      )
      
      coef_matrix <- fit_best$beta 
      selection_mat <- (coef_matrix != 0) * 1.0
      
      if(ncol(selection_mat) < n_lambda_seq) {
        selection_mat_padded <- matrix(0, nrow=p_vars, ncol=n_lambda_seq)
        selection_mat_padded[, 1:ncol(selection_mat)] <- selection_mat
        selection_mat <- selection_mat_padded
      }
      
      stability_counts <- stability_counts + selection_mat
      
    }, error = function(e) {
    })
  }
  
  if(message) close(pb)
  
  stability_probs <- stability_counts / B
  
  res <- list(
    stability_path = stability_probs,
    lambda = lambda
  )
  class(res) <- "StabSelect"
  return(res)
}



#' Stability Selection for Cox MDTL Enet Model
#'
#' @inheritParams cv.cox_MDTL_enet
#' @param B Integer. Number of bootstrap/subsampling replicates. Default is 50.
#' @param fraction_sample Numeric. Fraction of data to use for subsampling in each replicate. Default is 0.5.
#'
#' @return An object of class \code{StabSelect}.
#'
#' @examples
#' \dontrun{
#' data(ExampleData_highdim)
#' train_dat_highdim <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#'
#' eta_list <- generate_eta(method = "exponential", n = 10, max_eta = 10)
#'
#' mdtl.StabSelect <- cox_MDTL_enet.StabSelect(
#'   z = train_dat_highdim$z,
#'   delta = train_dat_highdim$status,
#'   time = train_dat_highdim$time,
#'   stratum = train_dat_highdim$stratum,
#'   beta = beta_external_highdim,
#'   vcov = NULL,
#'   etas = eta_list,
#'   alpha = 1,
#'   cv.criteria = "CIndex_pooled",
#'   B = 20,
#'   message = TRUE
#' )
#' plot(mdtl.StabSelect, threshold = 0.75)
#' }
#'
#' @export
cox_MDTL_enet.StabSelect <- function(z, delta, time, stratum = NULL,
                                     beta, vcov = NULL,
                                     etas = NULL, alpha = NULL,
                                     lambda = NULL, nlambda = 100, lambda.min.ratio = 0.1,
                                     nfolds = 5,
                                     cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
                                     c_index_stratum = NULL,
                                     message = FALSE, seed = NULL,
                                     B = 50,
                                     fraction_sample = 0.5,
                                     ...) {
  
  if (!is.null(seed)) set.seed(seed)
  cv.criteria <- match.arg(cv.criteria)
  
  z <- as.matrix(z)
  n_full <- nrow(z)
  p_vars <- ncol(z)
  var_names <- colnames(z)
  if(is.null(var_names)) var_names <- paste0("V", 1:p_vars)
  
  if (is.null(lambda.min.ratio)) {
    lambda.min.ratio <- ifelse(n_full < p_vars, 0.05, 1e-03)
  }
  
  if (is.null(lambda)) {
    if(message) message("Generating global lambda sequence using full data...")
    lambda_max_candidates <- numeric(length(etas))
    for(i in seq_along(etas)) {
      fit_tmp <- cox_MDTL_enet(z = z, delta = delta, time = time, stratum = stratum,
                               beta = beta, vcov = vcov, 
                               eta = etas[i], alpha = alpha,
                               nlambda = 5, lambda.min.ratio = lambda.min.ratio, ...)
      if(!is.null(fit_tmp$lambda)) lambda_max_candidates[i] <- max(fit_tmp$lambda)
    }
    global_lambda_max <- 1.1 * max(lambda_max_candidates)
    lambda <- exp(seq(log(global_lambda_max), log(global_lambda_max * lambda.min.ratio), length.out = nlambda))
  }
  
  n_lambda_seq <- length(lambda)
  stability_counts <- matrix(0, nrow = p_vars, ncol = n_lambda_seq)
  rownames(stability_counts) <- var_names
  
  if(message) {
    message(sprintf("Starting Stability Selection (MDTL) with B=%d...", B))
    pb <- txtProgressBar(min = 0, max = B, style = 3)
  }
  
  for (b in 1:B) {
    if(message) setTxtProgressBar(pb, b)
    
    sub_idx <- sort(sample(seq_len(n_full), size = floor(fraction_sample * n_full)))
    z_sub <- z[sub_idx, , drop = FALSE]
    delta_sub <- delta[sub_idx]
    time_sub <- time[sub_idx]
    stratum_sub <- if(!is.null(stratum)) stratum[sub_idx] else NULL
    
    tryCatch({
      cv_fit <- cv.cox_MDTL_enet(
        z = z_sub, delta = delta_sub, time = time_sub, stratum = stratum_sub,
        beta = beta, vcov = vcov,
        etas = etas, alpha = alpha,
        lambda = lambda,
        nfolds = nfolds,
        cv.criteria = cv.criteria,
        message = FALSE,
        ...
      )
      
      best_eta <- cv_fit$best$best_eta
      
      fit_best <- cox_MDTL_enet(
        z = z_sub, delta = delta_sub, time = time_sub, stratum = stratum_sub,
        beta = beta, vcov = vcov,
        eta = best_eta,
        alpha = alpha,
        lambda = lambda,
        ...
      )
      
      coef_matrix <- fit_best$beta
      selection_mat <- (coef_matrix != 0) * 1.0
      
      if(ncol(selection_mat) < n_lambda_seq) {
        selection_mat_padded <- matrix(0, nrow=p_vars, ncol=n_lambda_seq)
        selection_mat_padded[, 1:ncol(selection_mat)] <- selection_mat
        selection_mat <- selection_mat_padded
      }
      
      stability_counts <- stability_counts + selection_mat
      
    }, error = function(e) {
    })
  }
  
  if(message) close(pb)
  
  stability_probs <- stability_counts / B
  
  res <- list(
    stability_path = stability_probs,
    lambda = lambda
  )
  class(res) <- "StabSelect"
  return(res)
}