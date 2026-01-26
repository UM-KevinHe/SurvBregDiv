#' Plot Cross-Validation Results vs Eta (Biometrics-style)
#'
#' @description
#' Plots cross-validation performance across eta values for
#' \code{cv.coxkl}, \code{cv.coxkl_ridge}, or \code{cv.coxkl_enet} objects
#' in the Biometrics figure style. It displays the CV performance curve,
#' a baseline reference at eta=0, and marks the optimal eta.
#'
#' @param object A fitted cross-validation result of class \code{"cv.coxkl"},
#'   \code{"cv.coxkl_ridge"}, or \code{"cv.coxkl_enet"}.
#' @param line_color Color for the CV performance curve. Default is \code{"#7570B3"}.
#' @param baseline_color Color for the external baseline line. Default is \code{"#1B9E77"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object showing cross-validation performance versus \code{eta}.
#'
#' @examples
#' data(ExampleData_lowdim)
#' train_dat_lowdim <- ExampleData_lowdim$train
#' test_dat_lowdim <- ExampleData_lowdim$test
#' beta_external_lowdim <- ExampleData_lowdim$beta_external_good
#'
#' eta_list <- generate_eta(method = "exponential", n = 50, max_eta = 30)
#' cv.coxkl_est <- cv.coxkl(z = train_dat_lowdim$z,
#'                          delta = train_dat_lowdim$status,
#'                          time = train_dat_lowdim$time,
#'                          beta = beta_external_lowdim,
#'                          etas = eta_list,
#'                          criteria = "V&VH",
#'                          seed = 1)
#' cv.plot(cv.coxkl_est)
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid get_legend
#' @importFrom grid unit
#' @importFrom rlang .data
#' @export
cv.plot <- function(object,
                    line_color = "#7570B3",
                    baseline_color = "#1B9E77",
                    ...) {
  if (inherits(object, "cv.coxkl") || inherits(object, "cv.cox_MDTL")) {
    df <- object$internal_stat
    criteria <- object$criteria
  } else if (inherits(object, "cv.coxkl_ridge") ||
             inherits(object, "cv.cox_MDTL_ridge") ||
             inherits(object, "cv.coxkl_enet") ||
             inherits(object, "cv.cox_MDTL_enet")) {
    df <- object$integrated_stat.best_per_eta
    criteria <- object$criteria
  } else {
    stop("Object class not recognized for cv.plot.", call. = FALSE)
  }
  
  is_loss <- criteria %in% c("V&VH", "LinPred")
  ylab <- if (is_loss) "Loss" else "C Index"
  
  loss_candidates   <- c("VVH_Loss", "LinPred_Loss", "Loss", "loss")
  cindex_candidates <- c("CIndex_pooled", "CIndex_foldaverage", "CIndex", "cindex")
  candidates <- if (is_loss) loss_candidates else cindex_candidates
  
  metric_col <- NULL
  for (nm in candidates) if (nm %in% names(df)) { metric_col <- nm; break }
  if (is.null(metric_col)) {
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    num_cols <- setdiff(num_cols, c("eta", "lambda"))
    if (length(num_cols) == 0L) stop("Could not detect metric column in CV results.", call. = FALSE)
    metric_col <- num_cols[length(num_cols)]
  }
  
  df$eta <- as.numeric(df$eta)
  df <- df[order(df$eta), , drop = FALSE]
  df$metric <- as.numeric(df[[metric_col]])
  
  if (!any(df$eta == 0)) {
    idx0 <- which.min(abs(df$eta - 0))
  } else {
    idx0 <- which(df$eta == 0)[1]
  }
  baseline_val <- df$metric[idx0]
  baseline_eta <- df$eta[idx0]
  
  if (is_loss) {
    opt_idx <- which.min(df$metric)
  } else {
    opt_idx <- which.max(df$metric)
  }
  opt_eta <- df$eta[opt_idx]
  
  xmin <- min(df$eta, na.rm = TRUE)
  xmax <- max(df$eta, na.rm = TRUE)
  
  y_all <- c(df$metric, baseline_val)
  ylow  <- min(y_all, na.rm = TRUE) * 0.995
  yhigh <- max(y_all, na.rm = TRUE) * 1.005
  
  g_main <- ggplot(df, aes(x = .data$eta, y = .data$metric, group = 1)) +
    geom_line(linewidth = 1, color = line_color) +
    geom_point(size = 1.3, color = line_color) +
    geom_segment(
      data = data.frame(xmin = xmin, xmax = xmax, y = baseline_val),
      aes(x = .data$xmin, xend = .data$xmax, y = .data$y, yend = .data$y),
      inherit.aes = FALSE, color = baseline_color, linetype = "dotted", linewidth = 1
    ) +
    geom_point(
      data = data.frame(eta = baseline_eta, metric = baseline_val),
      aes(x = .data$eta, y = .data$metric),
      inherit.aes = FALSE, color = baseline_color, shape = 16, size = 2.4
    ) +
    geom_segment(
      data = data.frame(x = opt_eta),
      aes(x = .data$x, xend = .data$x, y = ylow, yend = yhigh),
      inherit.aes = FALSE, color = "#D95F02", linewidth = 1, linetype = "dashed"
    ) +
    labs(x = expression(eta), y = ylab) +
    coord_cartesian(ylim = c(ylow, yhigh)) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.length = unit(0.1, "cm"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(size = 14),
      legend.position = "none"
    )
  
  legend_df <- data.frame(
    x = rep(c(0, 1), 2),
    y = rep(1, 4),
    Method = factor(rep(c("Integrated", "Internal"), each = 2),
                    levels = c("Integrated", "Internal"))
  )
  
  g_legend <- ggplot(legend_df, aes(x = .data$x, y = .data$y, color = .data$Method, linetype = .data$Method)) +
    geom_line(linewidth = 1) +
    geom_point(
      data = subset(legend_df, Method == "Internal"),
      aes(x = 0.5, y = 1, color = Method),
      inherit.aes = FALSE, shape = 16, size = 2.4
    ) +
    scale_color_manual(values = c("Integrated" = line_color, "Internal" = baseline_color)) +
    scale_linetype_manual(values = c("Integrated" = "solid", "Internal" = "dotted")) +
    theme_void(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.2, "lines"),
      legend.key.height = unit(0.6, "lines")
    ) +
    guides(color = guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL),
           linetype = guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL))
  
  cowplot::plot_grid(cowplot::get_legend(g_legend), g_main, ncol = 1, rel_heights = c(0.08, 1))
}