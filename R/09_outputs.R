summarize_design_size_sensitivity <- function(sim) {
  require_output_packages()
  dplyr::as_tibble(sim) |>
    dplyr::group_by(experiment, design_id, n_seq, seq_len, n_total) |>
    dplyr::summarise(
      `Pearson iid` = mean(reject_naive, na.rm = TRUE),
      `RS first-order` = mean(reject_rs1, na.rm = TRUE),
      `RS second-order` = mean(reject_rs2, na.rm = TRUE),
      `Sequence Wald` = mean(reject_wald, na.rm = TRUE),
      `Mean design effect` = mean(lambda_mean, na.rm = TRUE),
      `Mean RS2 df` = mean(rs2_df, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()
}

plot_design_size_sensitivity <- function(sim) {
  require_output_packages()
  dd <- dplyr::as_tibble(sim) |>
    dplyr::select(experiment, design_id, n_seq, seq_len, n_total,
                  reject_naive, reject_rs2, reject_wald) |>
    tidyr::pivot_longer(
      cols = c(reject_naive, reject_rs2, reject_wald),
      names_to = "method_id", values_to = "reject"
    ) |>
    dplyr::mutate(
      method = dplyr::recode(
        method_id,
        reject_naive = "Pearson + iid reference",
        reject_rs2 = "second-order design correction",
        reject_wald = "sequence-aware Wald"
      )
    ) |>
    dplyr::group_by(experiment, design_id, n_seq, seq_len, n_total, method) |>
    dplyr::summarise(
      rate = mean(reject, na.rm = TRUE),
      R = dplyr::n(),
      se = sqrt(rate * (1 - rate) / R),
      .groups = "drop"
    )

  ggplot2::ggplot(dd, ggplot2::aes(x = seq_len, y = rate,
                                   linetype = method, shape = method)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, rate - 1.96 * se),
                   ymax = pmin(1, rate + 1.96 * se)),
      width = 8
    ) +
    ggplot2::scale_x_log10(breaks = sort(unique(dd$seq_len))) +
    ggplot2::facet_wrap(~experiment, scales = "free_x") +
    ggplot2::labs(
      x = "Sequence length L (log scale)",
      y = "Rejection rate at nominal 5%",
      linetype = NULL, shape = NULL,
      title = "The calibration problem is not an artifact of n = 100,000",
      subtitle = "Wrapped-Gaussian circular Markov construction with rho = 0.5"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

summarize_cvm_calibration <- function(sim) {
  require_output_packages()
  dplyr::as_tibble(sim) |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      `CvM + iid reference` = mean(reject_cvm_iid, na.rm = TRUE),
      `CvM + design correction` = mean(reject_cvm_rs2, na.rm = TRUE),
      `Mean CvM statistic` = mean(cvm_grid, na.rm = TRUE),
      `Mean RS2 df` = mean(cvm_rs2_df, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()
}

plot_cvm_calibration <- function(sim) {
  require_output_packages()
  dd <- dplyr::as_tibble(sim) |>
    dplyr::select(model, reject_cvm_iid, reject_cvm_rs2) |>
    tidyr::pivot_longer(-model, names_to = "method_id", values_to = "reject") |>
    dplyr::mutate(
      method = dplyr::recode(
        method_id,
        reject_cvm_iid = "CvM + iid reference",
        reject_cvm_rs2 = "CvM + design correction"
      )
    ) |>
    dplyr::group_by(model, method) |>
    dplyr::summarise(
      rate = mean(reject, na.rm = TRUE),
      R = dplyr::n(),
      se = sqrt(rate * (1 - rate) / R),
      .groups = "drop"
    )

  ggplot2::ggplot(dd, ggplot2::aes(x = method, y = rate)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, rate - 1.96 * se),
                   ymax = pmin(1, rate + 1.96 * se)),
      width = 0.1
    ) +
    ggplot2::facet_wrap(~model, ncol = 2) +
    ggplot2::labs(
      x = NULL,
      y = "Rejection rate at nominal 5%",
      title = "The same sampling issue appears for a continuous-significand CvM statistic",
      subtitle = "Second-order calibration uses the sequence-level covariance of the empirical CDF"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

summarize_marginal_power <- function(sim) {
  require_output_packages()
  dplyr::as_tibble(sim) |>
    dplyr::group_by(epsilon) |>
    dplyr::summarise(
      `Pearson iid` = mean(reject_pearson_iid, na.rm = TRUE),
      `Pearson RS2` = mean(reject_pearson_rs2, na.rm = TRUE),
      `Pearson Wald` = mean(reject_pearson_wald, na.rm = TRUE),
      `CvM iid` = mean(reject_cvm_iid, na.rm = TRUE),
      `CvM design` = mean(reject_cvm_rs2, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()
}

plot_marginal_power <- function(sim) {
  require_output_packages()
  dd <- dplyr::as_tibble(sim) |>
    dplyr::select(epsilon,
                  reject_pearson_iid, reject_pearson_rs2, reject_pearson_wald,
                  reject_cvm_iid, reject_cvm_rs2) |>
    tidyr::pivot_longer(-epsilon, names_to = "method_id", values_to = "reject") |>
    dplyr::mutate(
      family = ifelse(grepl("cvm", method_id), "Continuous-significand CvM", "First-digit Pearson"),
      method = dplyr::recode(
        method_id,
        reject_pearson_iid = "iid reference",
        reject_pearson_rs2 = "second-order design correction",
        reject_pearson_wald = "sequence-aware Wald",
        reject_cvm_iid = "iid reference",
        reject_cvm_rs2 = "second-order design correction"
      )
    ) |>
    dplyr::group_by(epsilon, family, method) |>
    dplyr::summarise(
      rate = mean(reject, na.rm = TRUE),
      R = dplyr::n(),
      se = sqrt(rate * (1 - rate) / R),
      .groups = "drop"
    )

  ggplot2::ggplot(dd, ggplot2::aes(x = epsilon, y = rate,
                                   linetype = method, shape = method)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, rate - 1.96 * se),
                   ymax = pmin(1, rate + 1.96 * se)),
      width = 0.0015
    ) +
    ggplot2::facet_wrap(~family, ncol = 2) +
    ggplot2::labs(
      x = expression(epsilon~"in"~f[epsilon](u)==1+epsilon*(2*u-1)),
      y = "Rejection probability",
      linetype = NULL, shape = NULL,
      title = "Design-aware calibration restores size without eliminating power",
      subtitle = "Smooth marginal alternatives are imposed on the same rho = 0.5 circular Markov dependence"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
