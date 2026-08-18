# Generators for the Benford equifinality benchmark (v2)
# All primary generators return a list with a matrix `u` of logarithmic
# significands in [0,1), with rows = sequences and columns = time points.
#
# v2 change: the PRIMARY latent-regime generator uses an exactly balanced
# number of Z=0 and Z=1 sequences. A separate random-composition generator is
# retained for the cluster-composition sensitivity experiment.

mod1 <- function(x) x %% 1

gen_iid_benford <- function(n_seq = 400L, seq_len = 250L) {
  u <- matrix(runif(n_seq * seq_len), nrow = n_seq, ncol = seq_len)
  list(model = "iid Benford", u = u, latent = NULL)
}

gen_markov_benford <- function(n_seq = 400L, seq_len = 250L, rho = 0.5) {
  stopifnot(rho > 0, rho < 1)
  sigma2 <- -log(rho) / (2 * pi^2)
  sigma <- sqrt(sigma2)

  u <- matrix(NA_real_, nrow = n_seq, ncol = seq_len)
  u[, 1L] <- runif(n_seq)
  if (seq_len >= 2L) {
    for (tt in 2:seq_len) {
      eps <- rnorm(n_seq, mean = 0, sd = sigma)
      u[, tt] <- mod1(u[, tt - 1L] + eps)
    }
  }

  list(
    model = "stochastic multiplicative",
    u = u,
    latent = NULL,
    rho = rho,
    sigma = sigma
  )
}

gen_rotation_benford <- function(n_seq = 400L, seq_len = 250L,
                                 theta = (sqrt(5) - 1) / 2) {
  u0 <- runif(n_seq)
  shifts <- (0:(seq_len - 1L)) * theta
  u <- outer(u0, shifts, "+") %% 1

  list(
    model = "deterministic rotation",
    u = u,
    latent = NULL,
    theta = theta
  )
}

# Helper used by both latent-regime designs.
.gen_latent_from_labels <- function(latent, seq_len, design) {
  n_seq <- length(latent)
  latent_obs <- rep(latent, each = seq_len)

  # Z=0 -> Beta(2,1), Z=1 -> Beta(1,2).
  shape1 <- ifelse(latent_obs == 0L, 2, 1)
  shape2 <- ifelse(latent_obs == 0L, 1, 2)
  vals <- rbeta(n_seq * seq_len, shape1 = shape1, shape2 = shape2)

  # rep(..., each=seq_len) is row-major; matrix() otherwise fills by columns.
  u <- matrix(vals, nrow = n_seq, ncol = seq_len, byrow = TRUE)

  list(
    model = "latent-regime mixture",
    u = u,
    latent = latent,
    mixture_design = design,
    pi_hat = mean(latent == 0L)
  )
}

# PRIMARY latent control: exactly half the sequences are in each regime.
# This removes random cluster-composition noise while retaining the key
# aggregation property: neither conditional regime is Benford, but the
# equal-weight pooled population is exactly uniform.
gen_latent_mixture_balanced <- function(n_seq = 400L, seq_len = 250L) {
  n_seq <- as.integer(n_seq)
  if (n_seq %% 2L != 0L) {
    stop("Balanced latent mixture requires an even n_seq; got ", n_seq, ".")
  }
  latent <- sample(rep(c(0L, 1L), each = n_seq / 2L), size = n_seq, replace = FALSE)
  .gen_latent_from_labels(latent, seq_len = seq_len, design = "balanced 50/50")
}

# Sensitivity control: sequence-level regimes are sampled independently.
# The SUPEROPOPULATION mixture is exactly 50/50 Benford, but a finite sample
# generally realizes pi_hat != 1/2, producing cluster-composition noise.
gen_latent_mixture_random <- function(n_seq = 400L, seq_len = 250L) {
  latent <- rbinom(n_seq, size = 1L, prob = 0.5)
  .gen_latent_from_labels(latent, seq_len = seq_len,
                          design = "random Bernoulli composition")
}

# Backward-compatible name: primary benchmark now means balanced composition.
gen_latent_mixture <- gen_latent_mixture_balanced

generate_model <- function(model, n_seq = 400L, seq_len = 250L,
                           rho = 0.5,
                           theta = (sqrt(5) - 1) / 2) {
  switch(
    model,
    iid = gen_iid_benford(n_seq, seq_len),
    markov = gen_markov_benford(n_seq, seq_len, rho = rho),
    rotation = gen_rotation_benford(n_seq, seq_len, theta = theta),
    mixture = gen_latent_mixture_balanced(n_seq, seq_len),
    stop("Unknown model: ", model)
  )
}
