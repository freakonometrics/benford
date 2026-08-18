# Design-aware calibration for first-digit Benford goodness-of-fit tests
# Pass 2: sequence-level covariance estimation, Rao-Scott-style moment
# corrections, and a sequence-aware Wald test.
BENFORD_DESIGN_VERSION <- "2026-08-17-pass2-design-aware-v1"

cluster_digit_proportions <- function(u) {
  u <- as.matrix(u); B <- nrow(u); L <- ncol(u)
  d <- matrix(first_digit_from_u(as.vector(u)), nrow = B, ncol = L)
  out <- matrix(0, nrow = B, ncol = 9L)
  for (j in 1:9) out[, j] <- rowMeans(d == j)
  colnames(out) <- paste0("d", 1:9); out
}

pearson_reduced_metric <- function(p = benford_first_digit_probs()) {
  p8 <- p[1:8]
  diag(1 / p8, nrow = 8L) + matrix(1 / p[9], nrow = 8L, ncol = 8L)
}
.symmetrize <- function(M) (M + t(M)) / 2
.symmetric_sqrt <- function(M, tol = 1e-12) {
  ee <- eigen(.symmetrize(M), symmetric = TRUE)
  cutoff <- tol * max(1, max(abs(ee$values))); vals <- pmax(ee$values, 0); vals[vals < cutoff] <- 0
  ee$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(ee$vectors)
}
.symmetric_pinv <- function(M, tol = 1e-10) {
  ee <- eigen(.symmetrize(M), symmetric = TRUE)
  cutoff <- tol * max(1, max(abs(ee$values))); keep <- ee$values > cutoff
  if (!any(keep)) return(list(inv = matrix(0, nrow(M), ncol(M)), rank = 0L, eigenvalues = ee$values))
  V <- ee$vectors[, keep, drop = FALSE]
  list(inv = V %*% diag(1 / ee$values[keep], nrow = sum(keep)) %*% t(V), rank = sum(keep), eigenvalues = ee$values)
}

estimate_sequence_omega <- function(obj, stratify_balanced = TRUE) {
  u <- as.matrix(obj$u); B <- nrow(u); L <- ncol(u)
  X <- cluster_digit_proportions(u)[, 1:8, drop = FALSE]
  use_strata <- isTRUE(stratify_balanced) && !is.null(obj$latent) && !is.null(obj$mixture_design) && identical(as.character(obj$mixture_design), "balanced 50/50")
  if (!use_strata) {
    if (B < 3L) stop("At least three independent sequences are required.")
    return(list(omega = L * stats::cov(X), p_hat8 = colMeans(X), n_clusters = B, n_strata = 1L, covariance_design = "sequence-cluster"))
  }
  z <- as.factor(obj$latent); lev <- levels(z); p_hat8 <- numeric(8L); omega <- matrix(0, 8L, 8L)
  for (h in lev) {
    ind <- which(z == h); if (length(ind) < 3L) stop("Each fixed stratum needs at least three sequences.")
    w <- length(ind) / B; Xh <- X[ind, , drop = FALSE]
    p_hat8 <- p_hat8 + w * colMeans(Xh); omega <- omega + L * w * stats::cov(Xh)
  }
  list(omega = omega, p_hat8 = p_hat8, n_clusters = B, n_strata = length(lev), covariance_design = "fixed-strata sequence-cluster")
}

design_effect_spectrum <- function(omega, p = benford_first_digit_probs(), tol = 1e-10) {
  Ahalf <- .symmetric_sqrt(pearson_reduced_metric(p)); M <- .symmetrize(Ahalf %*% omega %*% Ahalf)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values; ev[abs(ev) < tol] <- 0; pmax(ev, 0)
}

rao_scott_style <- function(pearson_stat, lambdas, alpha = 0.05, tol = 1e-10) {
  lam <- lambdas[lambdas > tol]; r <- length(lam)
  if (r == 0L) return(c(lambda_bar=NA,rs1_stat=NA,rs1_p=NA,rs2_scale=NA,rs2_df=NA,rs2_stat=NA,rs2_p=NA))
  lambda_bar <- mean(lam); rs1_stat <- pearson_stat/lambda_bar; rs1_p <- pchisq(rs1_stat,df=r,lower.tail=FALSE)
  scale2 <- sum(lam^2)/sum(lam); df2 <- sum(lam)^2/sum(lam^2); rs2_stat <- pearson_stat/scale2; rs2_p <- pchisq(rs2_stat,df=df2,lower.tail=FALSE)
  c(lambda_bar=lambda_bar,rs1_stat=rs1_stat,rs1_p=rs1_p,rs2_scale=scale2,rs2_df=df2,rs2_stat=rs2_stat,rs2_p=rs2_p)
}

sequence_wald_benford <- function(obj, p = benford_first_digit_probs(), stratify_balanced = TRUE, tol = 1e-10) {
  est <- estimate_sequence_omega(obj, stratify_balanced); n <- length(obj$u); diff <- est$p_hat8-p[1:8]; pinv <- .symmetric_pinv(est$omega,tol)
  stat <- as.numeric(n*crossprod(diff,pinv$inv%*%diff)); df <- pinv$rank; pv <- if(df>0) pchisq(stat,df=df,lower.tail=FALSE) else NA_real_
  list(statistic=stat,df=df,p_value=pv,reject_05=is.finite(pv)&&pv<.05,omega=est$omega,p_hat8=est$p_hat8,n_clusters=est$n_clusters,n_strata=est$n_strata,covariance_design=est$covariance_design)
}

design_aware_benford_tests <- function(obj, pearson_stat = NULL, p = benford_first_digit_probs(), stratify_balanced = TRUE, tol = 1e-10) {
  if (is.null(pearson_stat)) pearson_stat <- benford_chisq(obj$u)$statistic
  est <- estimate_sequence_omega(obj,stratify_balanced); lam <- design_effect_spectrum(est$omega,p,tol); rs <- rao_scott_style(pearson_stat,lam,tol=tol); w <- sequence_wald_benford(obj,p,stratify_balanced,tol)
  list(lambda=lam,lambda_mean=unname(rs["lambda_bar"]),lambda_max=max(lam),lambda_min_positive=if(any(lam>tol))min(lam[lam>tol])else NA_real_,rs1_stat=unname(rs["rs1_stat"]),rs1_p=unname(rs["rs1_p"]),reject_rs1_05=unname(rs["rs1_p"])<.05,rs2_scale=unname(rs["rs2_scale"]),rs2_df=unname(rs["rs2_df"]),rs2_stat=unname(rs["rs2_stat"]),rs2_p=unname(rs["rs2_p"]),reject_rs2_05=unname(rs["rs2_p"])<.05,wald_stat=w$statistic,wald_df=w$df,wald_p=w$p_value,reject_wald_05=w$reject_05,covariance_design=w$covariance_design,n_clusters=w$n_clusters,n_strata=w$n_strata)
}
