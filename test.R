# Quick smoke test for revision pass 3.
# Run from the reproducibility directory:
#   source("test.R")

source("R/01_generators.R")
source("R/02_statistics.R")
source("R/07_design_aware_tests.R")
source("R/03_simulation_engine_v5.R")
source("R/04_outputs.R")
source("R/05_theory.R")
source("R/06_multiplication_game.R")
source("R/08_extensions.R")
source("R/09_outputs.R")

cat("Engine: ", BENFORD_ENGINE_VERSION, "\n", sep = "")
cat("Design: ", BENFORD_DESIGN_VERSION, "\n", sep = "")

# 1. Inverse-CDF sanity check.
set.seed(1)
v <- runif(10000)
for (eps in c(0, 0.02, 0.1)) {
  u <- inv_cdf_linear_tilt(v, eps)
  err <- max(abs(cdf_linear_tilt(u, eps) - v))
  if (!is.finite(err) || err > 1e-10) stop("Tilt inverse-CDF check failed.")
}
cat("smooth marginal tilt inversion.\n")

# 2. Small design-size experiment.
grid_small <- data.frame(
  experiment = c("small A", "small B"),
  n_seq = c(20L, 40L),
  seq_len = c(20L, 10L),
  n_total = c(400L, 400L),
  design_id = c("B20_L20", "B40_L10"),
  stringsAsFactors = FALSE
)
ds <- run_design_size_sensitivity(
  R = 2L, design_grid = grid_small, rho = 0.5,
  seed = 101L, progress = FALSE
)
stopifnot(nrow(ds) == 4L)
stopifnot(all(c("reject_naive", "reject_rs2", "reject_wald") %in% names(ds)))
cat("B/L design sensitivity.\n")

# 3. CvM design-aware test on four small models.
cvm <- run_cvm_calibration_experiment(
  R = 2L, n_seq = 20L, seq_len = 20L,
  rho = 0.5, grid_n = 19L,
  seed = 202L, progress = FALSE
)
stopifnot(nrow(cvm) == 8L)
stopifnot(all(is.finite(cvm$cvm_grid)))
stopifnot(all(is.finite(cvm$cvm_rs2_df)))
cat("design-aware CvM calibration.\n")

# 4. Small power experiment.
pw <- run_marginal_power_experiment(
  epsilon_grid = c(0, 0.02),
  R = 2L, n_seq = 20L, seq_len = 20L,
  rho = 0.5, cvm_grid_n = 19L,
  seed = 303L, progress = FALSE
)
stopifnot(nrow(pw) == 4L)
stopifnot(all(c("reject_pearson_rs2", "reject_cvm_rs2") %in% names(pw)))
cat("level/power experiment.\n")

# 5. Positive control for the information boundary.
bd <- run_manipulation_boundary_demo(
  n_seq = 20L, seq_len = 20L,
  q_heaping = 0.05, seed = 404L
)
stopifnot(nrow(bd) == 3L)
cat("information-boundary positive control.\n\n")

cat("ALL SMOKE TESTS COMPLETED.\n")
