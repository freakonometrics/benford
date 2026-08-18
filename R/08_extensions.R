inv_cdf_linear_tilt <- function(v, epsilon = 0) {
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 || epsilon > 1) stop("epsilon must be a scalar in [0,1].")
  v <- pmin(1,pmax(0,v)); if(epsilon<1e-12)return(v); a<-1-epsilon; (-a+sqrt(a^2+4*epsilon*v))/(2*epsilon)
}
cdf_linear_tilt <- function(u, epsilon = 0) {if(length(epsilon)!=1L||!is.finite(epsilon)||epsilon<0||epsilon>1)stop("epsilon must be a scalar in [0,1].");u<-pmin(1,pmax(0,u));(1-epsilon)*u+epsilon*u^2}
pdf_linear_tilt <- function(u, epsilon = 0) {if(length(epsilon)!=1L||!is.finite(epsilon)||epsilon<0||epsilon>1)stop("epsilon must be a scalar in [0,1].");out<-1+epsilon*(2*u-1);out[u<0|u>1]<-0;out}
apply_linear_tilt <- function(obj,epsilon=0){out<-obj;out$u<-matrix(inv_cdf_linear_tilt(as.vector(obj$u),epsilon),nrow=nrow(obj$u),ncol=ncol(obj$u));out$epsilon<-epsilon;out$base_model<-obj$model;out$model<-paste0(obj$model," + marginal tilt");out}
linear_tilt_digit_probs <- function(epsilon=0){a<-log10(1:9);b<-log10(2:10);cdf_linear_tilt(b,epsilon)-cdf_linear_tilt(a,epsilon)}
apply_round_number_heaping <- function(obj,q=.02,width=.015,seed=NULL){if(!is.null(seed))set.seed(seed);if(q<0||q>1)stop("q must be in [0,1].");if(width<=0||width>=1)stop("width must be in (0,1).");out<-obj;u<-as.vector(obj$u);hit<-runif(length(u))<q;u[hit]<-runif(sum(hit),0,width);out$u<-matrix(u,nrow=nrow(obj$u),ncol=ncol(obj$u));out$heaping_q<-q;out$heaping_width<-width;out$model<-paste0(obj$model," + round-number heaping");out}

default_pass3_design_grid <- function(){a<-data.frame(experiment="Fixed B = 400; vary L and n",n_seq=400L,seq_len=c(10L,25L,50L,100L,250L),stringsAsFactors=FALSE);b<-data.frame(experiment="Fixed n = 100,000; vary B and L",n_seq=c(100L,200L,400L,1000L,2000L),seq_len=c(1000L,500L,250L,100L,50L),stringsAsFactors=FALSE);out<-rbind(a,b);out$n_total<-out$n_seq*out$seq_len;out$design_id<-sprintf("B%d_L%d",out$n_seq,out$seq_len);out}
run_design_size_sensitivity <- function(R=500L,design_grid=default_pass3_design_grid(),rho=.5,seed=20260901L,progress=interactive()){
  design_grid<-as.data.frame(design_grid);required<-c("experiment","n_seq","seq_len","n_total","design_id");if(!all(required%in%names(design_grid)))stop("design_grid must contain: ",paste(required,collapse=", "))
  out<-vector("list",nrow(design_grid)*R);idx<-1L
  for(jj in seq_len(nrow(design_grid))){B<-as.integer(design_grid$n_seq[jj]);L<-as.integer(design_grid$seq_len[jj]);if(B<10L)stop("Need at least 10 independent sequences for this experiment.");if(progress)message("Pass 3 design sensitivity: ",design_grid$design_id[jj]," (",jj,"/",nrow(design_grid),")")
    for(rr in seq_len(R)){set.seed(seed+100000L*jj+rr);obj<-gen_markov_benford(n_seq=B,seq_len=L,rho=rho);chi<-benford_chisq(obj$u);da<-design_aware_benford_tests(obj,pearson_stat=chi$statistic);out[[idx]]<-data.frame(rep=rr,experiment=design_grid$experiment[jj],design_id=design_grid$design_id[jj],n_seq=B,seq_len=L,n_total=B*L,rho=rho,chisq=chi$statistic,reject_naive=chi$p_value<.05,reject_rs1=isTRUE(da$reject_rs1_05),reject_rs2=isTRUE(da$reject_rs2_05),reject_wald=isTRUE(da$reject_wald_05),lambda_mean=da$lambda_mean,rs2_df=da$rs2_df,stringsAsFactors=FALSE);idx<-idx+1L}}
  do.call(rbind,out)}

cvm_default_grid <- function(m=49L){m<-as.integer(m);if(m<9L)stop("Use at least 9 interior grid points.");seq_len(m)/(m+1)}
sequence_cdf_matrix <- function(obj,grid=cvm_default_grid()){u<-as.matrix(obj$u);B<-nrow(u);X<-vapply(grid,function(g)rowMeans(u<=g),numeric(B));if(is.null(dim(X)))X<-matrix(X,nrow=B);colnames(X)<-sprintf("g%.5f",grid);X}
estimate_sequence_omega_cdf <- function(obj,grid=cvm_default_grid(),stratify_balanced=TRUE){
  X<-sequence_cdf_matrix(obj,grid);B<-nrow(X);L<-ncol(as.matrix(obj$u));fixed_strata<-isTRUE(stratify_balanced)&&!is.null(obj$latent)&&!is.null(obj$mixture_design)&&identical(as.character(obj$mixture_design),"balanced 50/50")
  if(!fixed_strata){if(B<3L)stop("At least three independent sequences are required.");return(list(omega=L*cov(X),F_hat=colMeans(X),n_clusters=B,covariance_design="sequence-cluster"))}
  z<-as.factor(obj$latent);lev<-levels(z);omega<-matrix(0,ncol(X),ncol(X));F_hat<-numeric(ncol(X));for(h in lev){ind<-which(z==h);if(length(ind)<3L)stop("Each fixed stratum needs at least three sequences.");w<-length(ind)/B;Xh<-X[ind,,drop=FALSE];F_hat<-F_hat+w*colMeans(Xh);omega<-omega+L*w*cov(Xh)};list(omega=omega,F_hat=F_hat,n_clusters=B,covariance_design="fixed-strata sequence-cluster")}
cvm_grid_statistic <- function(obj,grid=cvm_default_grid()){est<-sequence_cdf_matrix(obj,grid);F_hat<-colMeans(est);n<-length(obj$u);h<-1/(length(grid)+1);n*h*sum((F_hat-grid)^2)}
CVM_IID_CRITICAL_05 <- 0.46136

design_aware_cvm_test <- function(obj,grid=cvm_default_grid(),alpha=.05,stratify_balanced=TRUE,tol=1e-10){
  est<-estimate_sequence_omega_cdf(obj,grid,stratify_balanced);n<-length(obj$u);h<-1/(length(grid)+1);diff<-est$F_hat-grid;stat<-n*h*sum(diff^2);ev<-eigen(.symmetrize(est$omega),symmetric=TRUE,only.values=TRUE)$values;lam<-h*pmax(ev,0);lam<-lam[lam>tol]
  if(!length(lam))return(list(statistic=stat,naive_p_approx=NA_real_,reject_naive_05=stat>CVM_IID_CRITICAL_05,rs2_scale=NA_real_,rs2_df=NA_real_,rs2_p=NA_real_,reject_rs2_05=NA,lambda_sum=0,lambda_max=0,covariance_design=est$covariance_design))
  scale2<-sum(lam^2)/sum(lam);df2<-sum(lam)^2/sum(lam^2);rs2_stat<-stat/scale2;rs2_p<-pchisq(rs2_stat,df=df2,lower.tail=FALSE);naive_p<-pchisq(stat/(1/15),df=2.5,lower.tail=FALSE)
  list(statistic=stat,naive_p_approx=naive_p,reject_naive_05=stat>CVM_IID_CRITICAL_05,rs2_scale=scale2,rs2_df=df2,rs2_stat=rs2_stat,rs2_p=rs2_p,reject_rs2_05=rs2_p<alpha,lambda_sum=sum(lam),lambda_max=max(lam),covariance_design=est$covariance_design)}

run_cvm_calibration_experiment <- function(R=500L,n_seq=400L,seq_len=250L,rho=.5,theta=(sqrt(5)-1)/2,grid_n=49L,seed=20260902L,progress=interactive()){
  models<-c("iid","markov","rotation","mixture");grid<-cvm_default_grid(grid_n);out<-vector("list",R*length(models));idx<-1L
  for(rr in seq_len(R)){if(progress&&(rr==1L||rr%%25L==0L))message("Pass 3 CvM replication ",rr,"/",R);for(mm in seq_along(models)){set.seed(seed+10000L*rr+100L*mm);obj<-generate_model(models[mm],n_seq=n_seq,seq_len=seq_len,rho=rho,theta=theta);cv<-design_aware_cvm_test(obj,grid=grid);out[[idx]]<-data.frame(rep=rr,model_id=models[mm],model=obj$model,cvm_grid=cv$statistic,reject_cvm_iid=isTRUE(cv$reject_naive_05),cvm_rs2_scale=cv$rs2_scale,cvm_rs2_df=cv$rs2_df,cvm_rs2_p=cv$rs2_p,reject_cvm_rs2=isTRUE(cv$reject_rs2_05),cvm_lambda_sum=cv$lambda_sum,cvm_lambda_max=cv$lambda_max,stringsAsFactors=FALSE);idx<-idx+1L}}
  do.call(rbind,out)}

run_marginal_power_experiment <- function(epsilon_grid=c(0,.005,.01,.02,.03,.05),R=500L,n_seq=400L,seq_len=250L,rho=.5,cvm_grid_n=49L,seed=20260903L,progress=interactive()){
  epsilon_grid<-as.numeric(epsilon_grid);grid<-cvm_default_grid(cvm_grid_n);out<-vector("list",length(epsilon_grid)*R);idx<-1L
  for(ee in seq_along(epsilon_grid)){eps<-epsilon_grid[ee];if(progress)message("Pass 3 power: epsilon = ",eps);for(rr in seq_len(R)){set.seed(seed+100000L*ee+rr);base<-gen_markov_benford(n_seq=n_seq,seq_len=seq_len,rho=rho);obj<-apply_linear_tilt(base,eps);chi<-benford_chisq(obj$u);da<-design_aware_benford_tests(obj,pearson_stat=chi$statistic,stratify_balanced=FALSE);cv<-design_aware_cvm_test(obj,grid=grid,stratify_balanced=FALSE);out[[idx]]<-data.frame(rep=rr,epsilon=eps,n_seq=n_seq,seq_len=seq_len,n_total=n_seq*seq_len,rho=rho,pearson=chi$statistic,reject_pearson_iid=chi$p_value<.05,reject_pearson_rs2=isTRUE(da$reject_rs2_05),reject_pearson_wald=isTRUE(da$reject_wald_05),cvm_grid=cv$statistic,reject_cvm_iid=isTRUE(cv$reject_naive_05),reject_cvm_rs2=isTRUE(cv$reject_rs2_05),stringsAsFactors=FALSE);idx<-idx+1L}}
  do.call(rbind,out)}

run_manipulation_boundary_demo <- function(n_seq=400L,seq_len=250L,q_heaping=.02,heaping_width=.015,seed=20260904L){
  set.seed(seed);base<-gen_iid_benford(n_seq=n_seq,seq_len=seq_len);invisible<-base;invisible$model<-"integer-order magnitude shift";heaped<-apply_round_number_heaping(base,q=q_heaping,width=heaping_width,seed=seed+1L)
  summarize_one<-function(obj,label){chi<-benford_chisq(obj$u);data.frame(mechanism=label,pearson=chi$statistic,pearson_p=chi$p_value,cvm=cvm_uniform(obj$u),B_K=fourier_discrepancy(obj$u,K=5L),stringsAsFactors=FALSE)}
  rbind(summarize_one(base,"Benford reference"),summarize_one(invisible,"significand-preserving magnitude shift"),summarize_one(heaped,"round-number heaping"))}
