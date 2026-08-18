# Statistics for the Benford equifinality benchmark (v2)

benford_first_digit_probs <- function() {
  d <- 1:9
  log10(1 + 1 / d)
}

first_digit_from_u <- function(u) {
  # 10^u is in [1,10), so floor gives digits 1,...,9.
  pmin(9L, floor(10^u))
}

benford_chisq <- function(u) {
  d <- first_digit_from_u(as.vector(u))
  counts <- tabulate(d, nbins = 9L)
  p <- benford_first_digit_probs()
  n <- length(d)
  expected <- n * p
  stat <- sum((counts - expected)^2 / expected)
  p_value <- pchisq(stat, df = 8L, lower.tail = FALSE)

  list(
    statistic = unname(stat),
    p_value = unname(p_value),
    counts = counts,
    expected = expected
  )
}

cvm_uniform <- function(u) {
  x <- sort(as.vector(u))
  n <- length(x)
  i <- seq_len(n)
  1 / (12 * n) + sum((x - (2 * i - 1) / (2 * n))^2)
}

fourier_discrepancy <- function(u, K = 5L) {
  x <- as.vector(u)
  vals <- vapply(seq_len(K), function(k) {
    zbar <- mean(exp(2i * pi * k * x))
    Mod(zbar)^2
  }, numeric(1))
  sum(vals)
}

# Exact conditional expectation of the first circular transition moment after
# a random permutation within each sequence. For z_j=exp(2*pi*i*u_j), a
# random ordered pair of distinct observations has
# E[z_next * Conj(z_current)] = (|sum z_j|^2-L)/(L*(L-1)).
shuffle_transition_moment <- function(u) {
  L <- ncol(u)
  if (L < 2L) stop("Sequence length must be at least 2.")

  z <- exp(2i * pi * u)
  s <- rowSums(z)
  per_seq <- (Mod(s)^2 - L) / (L * (L - 1))
  mean(per_seq)
}

observed_transition_moment <- function(u) {
  L <- ncol(u)
  if (L < 2L) stop("Sequence length must be at least 2.")
  z <- exp(2i * pi * u)
  mean(z[, 2:L, drop = FALSE] * Conj(z[, 1:(L - 1L), drop = FALSE]))
}

circular_transition_stats <- function(u) {
  obs <- observed_transition_moment(u)
  null <- shuffle_transition_moment(u)
  centered <- obs - null

  c(
    C_raw = Mod(obs),
    C_shuffle = Mod(null),
    C_centered = Mod(centered),
    C_centered_real = Re(centered),
    C_centered_imag = Im(centered)
  )
}

summarize_benford_sample <- function(obj, K_fourier = 5L,
                                     include_design_aware = TRUE) {
  u <- obj$u
  chi <- benford_chisq(u)
  circ <- circular_transition_stats(u)

  out <- data.frame(
    model = obj$model,
    chisq = chi$statistic,
    p_chisq = chi$p_value,
    reject_05 = chi$p_value < 0.05,
    cvm = cvm_uniform(u),
    B_K = fourier_discrepancy(u, K = K_fourier),
    C_raw = unname(circ["C_raw"]),
    C_shuffle = unname(circ["C_shuffle"]),
    C_centered = unname(circ["C_centered"]),
    rho = if (!is.null(obj$rho)) obj$rho else NA_real_,
    theta = if (!is.null(obj$theta)) obj$theta else NA_real_,
    pi_hat = if (!is.null(obj$pi_hat)) obj$pi_hat else NA_real_,
    mixture_design = if (!is.null(obj$mixture_design)) obj$mixture_design else NA_character_,
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_design_aware)) {
    if (!exists("design_aware_benford_tests", mode = "function")) {
      stop("design_aware_benford_tests() is unavailable. Source R/07_design_aware_tests.R before running simulations.")
    }
    da <- design_aware_benford_tests(obj, pearson_stat = chi$statistic)
    out$lambda_mean <- da$lambda_mean
    out$lambda_max <- da$lambda_max
    out$lambda_min_positive <- da$lambda_min_positive
    out$rs1_stat <- da$rs1_stat
    out$rs1_p <- da$rs1_p
    out$reject_rs1_05 <- da$reject_rs1_05
    out$rs2_scale <- da$rs2_scale
    out$rs2_df <- da$rs2_df
    out$rs2_stat <- da$rs2_stat
    out$rs2_p <- da$rs2_p
    out$reject_rs2_05 <- da$reject_rs2_05
    out$wald_stat <- da$wald_stat
    out$wald_df <- da$wald_df
    out$wald_p <- da$wald_p
    out$reject_wald_05 <- da$reject_wald_05
    out$covariance_design <- da$covariance_design
  }

  out
}

# Rank-based AUC without an external dependency. By convention, AUC > 1/2
# means larger values tend to occur in y than in x.
auc_rank <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n0 <- length(x)
  n1 <- length(y)
  if (n0 == 0L || n1 == 0L) return(NA_real_)
  r <- rank(c(x, y), ties.method = "average")
  r1 <- sum(r[(n0 + 1L):(n0 + n1)])
  (r1 - n1 * (n1 + 1) / 2) / (n0 * n1)
}

signed_auc_effect <- function(x, y) 2 * auc_rank(x, y) - 1

bootstrap_auc_effect <- function(x, y, B = 1000L, seed = 1L) {
  set.seed(seed)
  est <- signed_auc_effect(x, y)
  boots <- replicate(B, {
    xb <- sample(x, length(x), replace = TRUE)
    yb <- sample(y, length(y), replace = TRUE)
    signed_auc_effect(xb, yb)
  })
  q <- quantile(boots, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(effect = est, lo = q[1], hi = q[2])
}
