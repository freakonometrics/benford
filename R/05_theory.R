
first_digit_intervals <- function() {
  data.frame(
    digit = 1:9,
    a = log10(1:9),
    b = log10(2:10),
    p = log10(1 + 1 / (1:9))
  )
}

latent_digit_probabilities <- function() {
  z <- first_digit_intervals()
  # Z=0: f_0(u)=2u; Z=1: f_1(u)=2(1-u)
  p0 <- z$b^2 - z$a^2
  p1 <- 2 * (z$b - z$a) - (z$b^2 - z$a^2)
  data.frame(digit = z$digit, p = z$p, p0 = p0, p1 = p1, q = p0 - p1)
}

expected_chisq_iid <- function() 8

# |\int_a^b exp(-2*pi*i*k*u)du|^2.
.interval_fourier_abs2 <- function(a, b, k) {
  (sin(pi * k * (b - a)) / (pi * k))^2
}

# Exact (up to a rapidly convergent Fourier truncation) expectation of the
# conventional first-digit Pearson statistic for the wrapped-Gaussian random
# walk, when each of B sequences starts from stationarity and has length L.
# The factor B cancels from E[T], so only L and rho matter.
expected_chisq_markov <- function(rho, seq_len = 250L, K_max = 100L) {
  stopifnot(rho >= 0, rho < 1, seq_len >= 2L)
  z <- first_digit_intervals()
  L <- as.integer(seq_len)
  ks <- seq_len(K_max)
  hs <- seq_len(L - 1L)

  ans <- 0
  for (j in seq_len(nrow(z))) {
    c2 <- .interval_fourier_abs2(z$a[j], z$b[j], ks)
    gamma <- vapply(hs, function(h) {
      # psi_k = rho^(k^2) for the wrapped Gaussian parametrization.
      2 * sum(c2 * rho^(ks^2 * h))
    }, numeric(1))

    var_one_sequence <-
      L * z$p[j] * (1 - z$p[j]) +
      2 * sum((L - hs) * gamma)

    ans <- ans + var_one_sequence / (L * z$p[j])
  }
  ans
}

# Length of overlap between A=[a,b) and A-shift on the unit circle.
.circular_interval_overlap <- function(a, b, shift) {
  shift <- shift %% 1
  lo <- (a - shift) %% 1
  hi <- (b - shift) %% 1

  overlap_linear <- function(l1, h1, l2, h2) {
    max(0, min(h1, h2) - max(l1, l2))
  }

  if (lo < hi) {
    overlap_linear(a, b, lo, hi)
  } else if (lo > hi) {
    overlap_linear(a, b, lo, 1) + overlap_linear(a, b, 0, hi)
  } else {
    0
  }
}

# Exact phase-averaged expectation for an irrational rotation with an
# independent uniform starting phase in each sequence.
expected_chisq_rotation <- function(theta = (sqrt(5) - 1) / 2,
                                    seq_len = 250L) {
  z <- first_digit_intervals()
  L <- as.integer(seq_len)
  hs <- seq_len(L - 1L)
  ans <- 0

  for (j in seq_len(nrow(z))) {
    cov_h <- vapply(hs, function(h) {
      ov <- .circular_interval_overlap(z$a[j], z$b[j], (h * theta) %% 1)
      ov - z$p[j]^2
    }, numeric(1))

    var_one_sequence <-
      L * z$p[j] * (1 - z$p[j]) +
      2 * sum((L - hs) * cov_h)

    # Small negative roundoff can occur when discrepancy is exceptionally low.
    var_one_sequence <- max(var_one_sequence, 0)
    ans <- ans + var_one_sequence / (L * z$p[j])
  }
  ans
}

# Exact conditional expectation of Pearson's statistic for a fixed realized
# fraction pi_hat of Z=0 sequences in the latent construction, assuming equal
# sequence lengths. At pi_hat=1/2 this reduces to the balanced benchmark.
expected_chisq_latent_composition <- function(pi_hat,
                                              n_total = 100000L) {
  pr <- latent_digit_probabilities()
  pi_hat <- as.numeric(pi_hat)
  delta <- pi_hat - 0.5

  variance_part <- sum(
    (pi_hat * pr$p0 * (1 - pr$p0) +
       (1 - pi_hat) * pr$p1 * (1 - pr$p1)) / pr$p
  )
  discrepancy_part <-
    as.numeric(n_total) * delta^2 * sum(pr$q^2 / pr$p)

  variance_part + discrepancy_part
}

expected_chisq_latent_balanced <- function() {
  expected_chisq_latent_composition(0.5, n_total = 1L)
}

latent_imbalance_coefficient <- function() {
  pr <- latent_digit_probabilities()
  sum(pr$q^2 / pr$p)
}

make_chisq_theory_table <- function(sim, rho, theta, seq_len) {
  mc <- tapply(sim$chisq, sim$model_id, mean)
  theory <- c(
    iid = expected_chisq_iid(),
    markov = expected_chisq_markov(rho = rho, seq_len = seq_len),
    rotation = expected_chisq_rotation(theta = theta, seq_len = seq_len),
    mixture = expected_chisq_latent_balanced()
  )
  labels <- c(
    iid = "iid Benford",
    markov = "stochastic multiplicative",
    rotation = "deterministic rotation",
    mixture = "latent-regime mixture"
  )
  ids <- names(theory)
  data.frame(
    Construction = unname(labels[ids]),
    `Theoretical E[chi-square]` = as.numeric(theory),
    `Monte Carlo mean` = as.numeric(mc[ids]),
    `MC minus theory` = as.numeric(mc[ids] - theory),
    check.names = FALSE
  )
}
