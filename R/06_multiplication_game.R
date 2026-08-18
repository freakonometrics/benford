
.game_mod1 <- function(x) x %% 1

draw_game_strategy <- function(n, strategy) {
  switch(strategy,
    haar = runif(n), beta_low = rbeta(n,2,5), beta_high = rbeta(n,5,2),
    wrapped_normal = .game_mod1(rnorm(n,mean=.22,sd=.07)),
    two_point = sample(c(.08,.58),n,replace=TRUE), stop("Unknown game strategy: ",strategy))
}

game_strategy_labels <- c(haar="Haar / Benford",beta_low="Beta(2,5) on log scale",beta_high="Beta(5,2) on log scale",wrapped_normal="concentrated wrapped normal",two_point="two-point strategy")
game_win_upper_default <- log10(4)
game_outcome <- function(a,b,win_upper=game_win_upper_default) (.game_mod1(a+b)<win_upper)

run_static_multiplication_game <- function(n_per_cell=50000L,win_upper=game_win_upper_default,seed=20260831L) {
  set.seed(seed); ids <- names(game_strategy_labels); out <- vector("list",length(ids)^2); z <- 1L
  for(r in ids) for(c in ids) {
    a <- draw_game_strategy(n_per_cell,r); b <- draw_game_strategy(n_per_cell,c); y <- game_outcome(a,b,win_upper); phat <- mean(y); se <- sqrt(phat*(1-phat)/n_per_cell)
    out[[z]] <- data.frame(row_strategy=r,column_strategy=c,row_label=unname(game_strategy_labels[r]),column_label=unname(game_strategy_labels[c]),win_probability=phat,se=se,target=win_upper,deviation=phat-win_upper,stringsAsFactors=FALSE); z <- z+1L
  }
  do.call(rbind,out)
}

summarize_haar_game_rows <- function(game_matrix) {
  a <- game_matrix[game_matrix$row_strategy=="haar",]; a$case <- "Player 1 uses Haar"; a$other_strategy <- a$column_label
  b <- game_matrix[game_matrix$column_strategy=="haar",]; b$case <- "Player 2 uses Haar"; b$other_strategy <- b$row_label
  out <- rbind(a,b); out[,c("case","other_strategy","win_probability","se","target","deviation")]
}

run_strategic_washout <- function(n=100000L,seed=20260832L) {
  set.seed(seed); u0 <- rbeta(n,5,2); b0 <- ifelse(u0<.5,.13,.67); a0 <- runif(n); u1 <- .game_mod1(u0+a0+b0)
  data <- rbind(data.frame(u=u0,stage="initial state"),data.frame(u=u1,stage="after one Haar move"))
  stat_one <- function(u){chi<-benford_chisq(u);data.frame(chisq=chi$statistic,p_chisq=chi$p_value,cvm=cvm_uniform(u),first_harmonic=Mod(mean(exp(2i*pi*u))))}
  tab <- rbind(stat_one(u0),stat_one(u1)); tab$stage <- c("initial state","after one Haar move"); tab <- tab[,c("stage","chisq","p_chisq","cvm","first_harmonic")]
  list(data=data,summary=tab)
}

run_nonhaar_repeated_game <- function(n_paths=50000L,T=25L,mu_a=.12,sigma_a=.08,mu_b=.47,sigma_b=.08,u0=.03,seed=20260833L) {
  set.seed(seed); u <- rep(u0%%1,n_paths); eta <- exp(-2*pi^2*(sigma_a^2+sigma_b^2)); rows <- vector("list",T+1L)
  record <- function(t,u)data.frame(t=t,first_harmonic=Mod(mean(exp(2i*pi*u))),theory=eta^t,stringsAsFactors=FALSE)
  rows[[1L]] <- record(0L,u)
  if(T>=1L) for(tt in seq_len(T)){a<-.game_mod1(rnorm(n_paths,mu_a,sigma_a));b<-.game_mod1(rnorm(n_paths,mu_b,sigma_b));u<-.game_mod1(u+a+b);rows[[tt+1L]]<-record(tt,u)}
  list(trajectory=do.call(rbind,rows),parameters=data.frame(parameter=c("mu_A","sigma_A","mu_B","sigma_B","eta"),value=c(mu_a,sigma_a,mu_b,sigma_b,eta),stringsAsFactors=FALSE))
}

plot_static_game_matrix <- function(game_matrix) {
  require_output_packages(); ggplot2::ggplot(game_matrix,ggplot2::aes(x=column_label,y=row_label,fill=deviation))+ggplot2::geom_tile()+ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",win_probability)),size=3)+ggplot2::labs(x="Player 2 strategy",y="Player 1 strategy",fill="Deviation\nfrom value",title="Haar randomization pins the multiplication-game payoff",subtitle="Cells report P(product significand has first digit 1, 2, or 3); game value = log10(4)")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30,hjust=1))
}
plot_strategic_washout <- function(washout) {
  require_output_packages(); grid<-data.frame(u=c(0,1),density=c(1,1)); ggplot2::ggplot(washout$data,ggplot2::aes(x=u))+ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)),bins=60,alpha=.28)+ggplot2::geom_line(data=grid,ggplot2::aes(x=u,y=density),inherit.aes=FALSE,linewidth=.8)+ggplot2::facet_wrap(~stage,ncol=2)+ggplot2::labs(x=expression(U=="{"*log[10](X)*"}"),y="Density",title="One Haar move washes out an arbitrary current-state response",subtitle="The opponent uses a deterministic state-dependent action; the solid line is Uniform(0,1)")+ggplot2::theme_minimal(base_size=11)
}
plot_nonhaar_game_mixing <- function(obj) {
  require_output_packages(); dd<-obj$trajectory; long<-rbind(data.frame(t=dd$t,value=dd$first_harmonic,series="simulation"),data.frame(t=dd$t,value=dd$theory,series="theory")); ggplot2::ggplot(long,ggplot2::aes(x=t,y=value,linetype=series))+ggplot2::geom_line(linewidth=.8)+ggplot2::geom_point(data=long[long$series=="simulation",],size=1.2)+ggplot2::scale_y_log10()+ggplot2::labs(x="Round t",y="Magnitude of first circular Fourier coefficient",linetype=NULL,title="Repeated non-Haar play can mix toward the Benford invariant law",subtitle="Both players use concentrated wrapped-normal log-actions; neither action distribution is Haar")+ggplot2::theme_minimal(base_size=11)
}
