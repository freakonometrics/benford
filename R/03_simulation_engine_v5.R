# Monte Carlo engine for the Benford equifinality benchmark (pass 2)
BENFORD_ENGINE_VERSION <- "2026-08-17-v5.1-design-aware-robust"

model_labels <- c(
  iid = "iid Benford",
  markov = "stochastic multiplicative",
  rotation = "deterministic rotation",
  mixture = "latent-regime mixture"
)

# Robust pass-2 wrapper.  The design-aware quantities are computed here even
# if an older copy of 02_statistics.R is accidentally sourced.  This avoids
# coupling the Monte Carlo engine to a particular version of
# summarize_benford_sample().
.ensure_design_aware_summary <- function(s, obj) {
  required <- c(
    "lambda_mean", "lambda_max", "lambda_min_positive",
    "rs1_stat", "rs1_p", "reject_rs1_05",
    "rs2_scale", "rs2_df", "rs2_stat", "rs2_p", "reject_rs2_05",
    "wald_stat", "wald_df", "wald_p", "reject_wald_05",
    "covariance_design"
  )

  if (all(required %in% names(s))) return(s)

  if (!exists("design_aware_benford_tests", mode = "function")) {
    stop(
      "Pass-2 design-aware functions are unavailable. ",
      "Source R/07_design_aware_tests.R before R/03_simulation_engine_v5.R."
    )
  }

  pearson_stat <- if ("chisq" %in% names(s) && length(s$chisq)) {
    as.numeric(s$chisq[1L])
  } else {
    benford_chisq(obj$u)$statistic
  }

  da <- design_aware_benford_tests(obj, pearson_stat = pearson_stat)
  for (nm in required) s[[nm]] <- da[[nm]]
  s
}

# Safe scalar extractor.  It turns a missing/NULL field into a typed NA rather
# than producing the opaque "replacement has length zero" error.
.scalar1 <- function(x, name, type = c("numeric", "logical", "character")) {
  type <- match.arg(type)
  z <- x[[name]]
  if (is.null(z) || length(z) == 0L) {
    return(switch(type, numeric = NA_real_, logical = NA, character = NA_character_))
  }
  z <- z[1L]
  switch(type, numeric = as.numeric(z), logical = as.logical(z), character = as.character(z))
}

.empty_primary_results <- function(n) {
  data.frame(
    rep = rep(NA_integer_, n),
    model_id = rep(NA_character_, n),
    model = rep(NA_character_, n),
    chisq = rep(NA_real_, n),
    p_chisq = rep(NA_real_, n),
    reject_05 = rep(NA, n),
    cvm = rep(NA_real_, n),
    B_K = rep(NA_real_, n),
    C_raw = rep(NA_real_, n),
    C_shuffle = rep(NA_real_, n),
    C_centered = rep(NA_real_, n),
    rho = rep(NA_real_, n),
    theta = rep(NA_real_, n),
    pi_hat = rep(NA_real_, n),
    mixture_design = rep(NA_character_, n),
    lambda_mean = rep(NA_real_, n),
    lambda_max = rep(NA_real_, n),
    lambda_min_positive = rep(NA_real_, n),
    rs1_stat = rep(NA_real_, n),
    rs1_p = rep(NA_real_, n),
    reject_rs1_05 = rep(NA, n),
    rs2_scale = rep(NA_real_, n),
    rs2_df = rep(NA_real_, n),
    rs2_stat = rep(NA_real_, n),
    rs2_p = rep(NA_real_, n),
    reject_rs2_05 = rep(NA, n),
    wald_stat = rep(NA_real_, n),
    wald_df = rep(NA_real_, n),
    wald_p = rep(NA_real_, n),
    reject_wald_05 = rep(NA, n),
    covariance_design = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

.write_summary_row <- function(out, i, s, rep_id, model_id) {
  stopifnot(is.data.frame(s), nrow(s) == 1L)
  out$rep[i] <- as.integer(rep_id)
  out$model_id[i] <- as.character(model_id)
  out$model[i] <- .scalar1(s, "model", "character")
  out$chisq[i] <- .scalar1(s, "chisq")
  out$p_chisq[i] <- .scalar1(s, "p_chisq")
  out$reject_05[i] <- .scalar1(s, "reject_05", "logical")
  out$cvm[i] <- .scalar1(s, "cvm")
  out$B_K[i] <- .scalar1(s, "B_K")
  out$C_raw[i] <- .scalar1(s, "C_raw")
  out$C_shuffle[i] <- .scalar1(s, "C_shuffle")
  out$C_centered[i] <- .scalar1(s, "C_centered")
  out$rho[i] <- .scalar1(s, "rho")
  out$theta[i] <- .scalar1(s, "theta")
  out$pi_hat[i] <- .scalar1(s, "pi_hat")
  out$mixture_design[i] <- .scalar1(s, "mixture_design", "character")
  out$lambda_mean[i] <- .scalar1(s, "lambda_mean")
  out$lambda_max[i] <- .scalar1(s, "lambda_max")
  out$lambda_min_positive[i] <- .scalar1(s, "lambda_min_positive")
  out$rs1_stat[i] <- .scalar1(s, "rs1_stat")
  out$rs1_p[i] <- .scalar1(s, "rs1_p")
  out$reject_rs1_05[i] <- .scalar1(s, "reject_rs1_05", "logical")
  out$rs2_scale[i] <- .scalar1(s, "rs2_scale")
  out$rs2_df[i] <- .scalar1(s, "rs2_df")
  out$rs2_stat[i] <- .scalar1(s, "rs2_stat")
  out$rs2_p[i] <- .scalar1(s, "rs2_p")
  out$reject_rs2_05[i] <- .scalar1(s, "reject_rs2_05", "logical")
  out$wald_stat[i] <- .scalar1(s, "wald_stat")
  out$wald_df[i] <- .scalar1(s, "wald_df")
  out$wald_p[i] <- .scalar1(s, "wald_p")
  out$reject_wald_05[i] <- .scalar1(s, "reject_wald_05", "logical")
  out$covariance_design[i] <- .scalar1(s, "covariance_design", "character")
  out
}


run_primary_simulations <- function(R = 200L,
                                    n_seq = 400L,
                                    seq_len = 250L,
                                    rho = 0.5,
                                    theta = (sqrt(5) - 1) / 2,
                                    K_fourier = 5L,
                                    seed = 20260817L,
                                    progress = interactive()) {
  if (n_seq %% 2L != 0L) stop("Primary v2 benchmark requires even n_seq.")
  models <- names(model_labels)
  out <- .empty_primary_results(as.integer(R) * length(models))
  idx <- 1L

  for (rr in seq_len(R)) {
    if (progress && (rr == 1L || rr %% 10L == 0L)) {
      message("Primary v2 replication ", rr, "/", R)
    }

    for (mm in seq_along(models)) {
      set.seed(seed + 10000L * rr + 100L * mm)
      obj <- generate_model(
        models[mm], n_seq = n_seq, seq_len = seq_len,
        rho = rho, theta = theta
      )
      s <- summarize_benford_sample(obj, K_fourier = K_fourier)
      s <- .ensure_design_aware_summary(s, obj)
      out <- .write_summary_row(out, idx, s, rr, models[mm])
      idx <- idx + 1L
    }
  }

  out$model_id <- factor(out$model_id, levels = models)
  out$model <- factor(out$model, levels = unname(model_labels))
  out
}

run_rho_sensitivity <- function(rho_grid = c(0, 0.1, 0.25, 0.5, 0.75, 0.9),
                                R = 100L,
                                n_seq = 400L,
                                seq_len = 250L,
                                K_fourier = 5L,
                                seed = 20260818L,
                                progress = interactive()) {
  out <- .empty_primary_results(length(rho_grid) * as.integer(R))
  idx <- 1L

  for (jj in seq_along(rho_grid)) {
    rho <- rho_grid[jj]
    for (rr in seq_len(R)) {
      if (progress && rr == 1L) message("rho = ", rho)
      set.seed(seed + 100000L * jj + rr)

      if (rho == 0) {
        obj <- gen_iid_benford(n_seq = n_seq, seq_len = seq_len)
        obj$model <- "stochastic multiplicative"
      } else {
        obj <- gen_markov_benford(n_seq = n_seq, seq_len = seq_len, rho = rho)
      }

      s <- summarize_benford_sample(obj, K_fourier = K_fourier)
      s <- .ensure_design_aware_summary(s, obj)
      s$rho <- rho
      out <- .write_summary_row(out, idx, s, rr, "markov")
      out$rho[idx] <- rho
      idx <- idx + 1L
    }
  }

  out$model_id <- factor(out$model_id, levels = "markov")
  out$model <- factor(out$model, levels = "stochastic multiplicative")
  out
}

# Separate experiment isolating finite-sample cluster-composition noise.
# Both designs use the same two non-Benford conditional laws. The balanced
# design fixes pi_hat=1/2; the random design draws sequence-level Z iid and
# therefore generally has pi_hat != 1/2 in a realized dataset.
run_mixture_composition_simulations <- function(R = 200L,
                                                n_seq = 400L,
                                                seq_len = 250L,
                                                K_fourier = 5L,
                                                seed = 20260822L,
                                                progress = interactive()) {
  if (n_seq %% 2L != 0L) stop("Balanced composition control requires even n_seq.")

  n_rows <- 2L * as.integer(R)
  out <- data.frame(
    rep = rep(NA_integer_, n_rows),
    design = rep(NA_character_, n_rows),
    pi_hat = rep(NA_real_, n_rows),
    imbalance = rep(NA_real_, n_rows),
    chisq = rep(NA_real_, n_rows),
    p_chisq = rep(NA_real_, n_rows),
    reject_05 = rep(NA, n_rows),
    cvm = rep(NA_real_, n_rows),
    B_K = rep(NA_real_, n_rows),
    C_centered = rep(NA_real_, n_rows),
    lambda_mean = rep(NA_real_, n_rows),
    rs1_p = rep(NA_real_, n_rows),
    reject_rs1_05 = rep(NA, n_rows),
    rs2_scale = rep(NA_real_, n_rows),
    rs2_df = rep(NA_real_, n_rows),
    rs2_p = rep(NA_real_, n_rows),
    reject_rs2_05 = rep(NA, n_rows),
    wald_p = rep(NA_real_, n_rows),
    reject_wald_05 = rep(NA, n_rows),
    covariance_design = rep(NA_character_, n_rows),
    stringsAsFactors = FALSE
  )

  idx <- 1L
  for (rr in seq_len(R)) {
    if (progress && (rr == 1L || rr %% 20L == 0L)) {
      message("Mixture-composition replication ", rr, "/", R)
    }

    generators <- list(
      "balanced 50/50" = gen_latent_mixture_balanced,
      "random Bernoulli composition" = gen_latent_mixture_random
    )

    for (jj in seq_along(generators)) {
      set.seed(seed + 10000L * rr + 100L * jj)
      obj <- generators[[jj]](n_seq = n_seq, seq_len = seq_len)
      s <- summarize_benford_sample(obj, K_fourier = K_fourier)
      s <- .ensure_design_aware_summary(s, obj)

      out$rep[idx] <- rr
      out$design[idx] <- names(generators)[jj]
      out$pi_hat[idx] <- obj$pi_hat
      out$imbalance[idx] <- abs(obj$pi_hat - 0.5)
      out$chisq[idx] <- s$chisq
      out$p_chisq[idx] <- s$p_chisq
      out$reject_05[idx] <- s$reject_05
      out$cvm[idx] <- s$cvm
      out$B_K[idx] <- s$B_K
      out$C_centered[idx] <- s$C_centered
      out$lambda_mean[idx] <- s$lambda_mean
      out$rs1_p[idx] <- s$rs1_p
      out$reject_rs1_05[idx] <- s$reject_rs1_05
      out$rs2_scale[idx] <- s$rs2_scale
      out$rs2_df[idx] <- s$rs2_df
      out$rs2_p[idx] <- s$rs2_p
      out$reject_rs2_05[idx] <- s$reject_rs2_05
      out$wald_p[idx] <- s$wald_p
      out$reject_wald_05[idx] <- s$reject_wald_05
      out$covariance_design[idx] <- s$covariance_design
      idx <- idx + 1L
    }
  }

  out$design <- factor(
    out$design,
    levels = c("balanced 50/50", "random Bernoulli composition")
  )
  out
}

make_representative_samples <- function(n_seq = 400L,
                                        seq_len = 250L,
                                        rho = 0.5,
                                        theta = (sqrt(5) - 1) / 2,
                                        seed = 20260819L) {
  if (n_seq %% 2L != 0L) stop("Representative v2 benchmark requires even n_seq.")
  models <- names(model_labels)
  out <- vector("list", length(models))
  names(out) <- models
  for (mm in seq_along(models)) {
    set.seed(seed + 100L * mm)
    out[[mm]] <- generate_model(
      models[mm], n_seq = n_seq, seq_len = seq_len,
      rho = rho, theta = theta
    )
  }
  out
}

run_manipulation_experiment <- function(n = 100000L,
                                        q = 0.10,
                                        shift_orders = 2L,
                                        K_values = 0:4,
                                        seed = 20260820L,
                                        K_fourier = 5L) {
  set.seed(seed)
  u <- runif(n)
  K <- sample(K_values, n, replace = TRUE)
  manipulated <- rbinom(n, size = 1L, prob = q)

  log10_x <- K + u
  log10_x_manip <- K + shift_orders * manipulated + u
  u_manip <- u

  before <- matrix(u, nrow = 1L)
  after <- matrix(u_manip, nrow = 1L)

  list(
    data = data.frame(
      u = u,
      K = K,
      manipulated = manipulated,
      log10_x = log10_x,
      log10_x_manip = log10_x_manip
    ),
    before = c(
      chisq = benford_chisq(before)$statistic,
      cvm = cvm_uniform(before),
      B_K = fourier_discrepancy(before, K = K_fourier)
    ),
    after = c(
      chisq = benford_chisq(after)$statistic,
      cvm = cvm_uniform(after),
      B_K = fourier_discrepancy(after, K = K_fourier)
    ),
    exact_u_identity = identical(u, u_manip)
  )
}
