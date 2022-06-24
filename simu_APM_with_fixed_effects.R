## janehoe@mail.ustc.edu.cn
## 06/25/2022

## Function to run a great deal of simulations for APM model on case-control studies 
## due to additive, parent-of-origin and maternal genetic effects with fixed effects.

## The functions follows the work of Laurin et al. (2018) (https://github.com/amatrhr/g-remladp)

Cc_generate <- function(dir,times,alpha0= 0.25,alpha1=1,alpha2=0.5,p_covar=0.5,var_x=1,hsqA=0.4,hsqP=0.05,hsqM=0.1,K=0.1,P=0.5,nsnps=10000,nobs){
  # dir: the directory of files
  # times: simulation id
  # hsqA: the true heritability of additive effects
  # hsqP: the true heritability of parent-of-origin effects
  # hsqM: the true heritability of maternal genetic effects
  # K: the prevalence
  # P: the case proportion in the study
  # nsnps: the number of SNPs
  # nobs: the number of children

  simname <- paste0("bin_",K)
  
  set.seed(times)
  aSD <-  function(pv, qv){
    sqrt(2 * pv * qv)
  }
  
  Dmean <- function(pv, qv){
    2 * pv^2
  }
  
  DSD <- function(pv, qv){
    sqrt(4 * pv^2 * qv^2)
  }
  
  # generate the MAFs and genetic effects
  true_EAFs <- runif(nsnps,0.1,0.5)
  betaA <- rnorm(n =  nsnps, mean = 0, sd = sqrt(hsqA / nsnps))
  betaP <- rnorm(n =  nsnps, mean = 0, sd = sqrt(hsqP / nsnps))
  betaM <- rnorm(n =  nsnps, mean = 0, sd = sqrt(hsqM / nsnps))
  
  write.table(cbind(betaA,betaP,betaM),paste0(dir,simname,"_APM_effects_",times,".txt"),col.names = c("A","P","M"),row.names = F)
 
  ############### sample selection ###########################
  case_mat_ped_snps <- matrix(0,ncol=nsnps)
  case_pat_ped_snps <- matrix(0,ncol=nsnps)
  case_mat_ped2_snps <- matrix(0,ncol=nsnps)
  control_mat_ped_snps <- matrix(0,ncol=nsnps)
  control_pat_ped_snps <- matrix(0,ncol=nsnps)
  control_mat_ped2_snps <- matrix(0,ncol=nsnps)
  case_covar_bin <- case_covar_conti <- control_covar_bin <- control_covar_conti <- c()
  
  
  tag <- 1
  while(tag <= nobs ){
    
    mat_ped_snps <- sapply(true_EAFs, rbinom, n = 1, size = 1) # the maternal-inherited SNP of the child
    pat_ped_snps <- sapply(true_EAFs, rbinom, n = 1, size = 1) # the paternal-inherited SNP of the child
    mat_ped2_snps <- sapply(true_EAFs, rbinom, n = 1, size = 1) # the other SNP of the mother which didn't pass to her child
    
    simulate_one_struct <- list()
    
    simulate_one_struct$mat_ped_snps <- mat_ped_snps
    simulate_one_struct$pat_ped_snps <- pat_ped_snps
    simulate_one_struct$mat_ped2_snps <- mat_ped2_snps
    
    p_m <- p <- true_EAFs
    q_m <- q <- 1 - p
    
    ## generate the genotypes of mother-child duos and recode them due to ADD, POE and MGE
    Aeps <- (simulate_one_struct$mat_ped_snps + simulate_one_struct$pat_ped_snps)
    Aeps <- matrix(Aeps,nrow=1)
    Peps <- (simulate_one_struct$mat_ped_snps - simulate_one_struct$pat_ped_snps)
    Peps <- matrix(Peps,nrow=1)
    Msnp <- (simulate_one_struct$mat_ped_snps + simulate_one_struct$mat_ped2_snps)
    Msnp <- matrix(Msnp,nrow=1)
    Meps <- Msnp * matrix(1, nrow = nrow(Msnp)) %*% matrix(2*p_m, nrow = 1) - 2*(Msnp > 1)
    
    ## standardise the codings
    Aeps <- sweep( sweep(Aeps, 2, 2 * p), 2, aSD(p,q), FUN = "/") 
    Peps <- sweep( Peps, 2, aSD(p,q), FUN = "/") 
    Meps <- sweep( sweep(Meps, 2, Dmean(p_m,q_m)), 2, DSD(p_m,q_m), FUN = "/") 
    
    ## generate the phenotypes
    thres <- qnorm(1-K,sd=sqrt(1+p_covar*(1-p_covar)*alpha1^2+var_x*alpha2^2))  ## threshold for the liability
    covar_bin <-  rbinom(1,size = 1,prob = p_covar)
    covar_conti <- rnorm(1,sd=sqrt(var_x))
    thres_single <- thres-(alpha1 * I(covar_bin==1)+alpha2 * covar_conti)
    
    liability <- alpha0 + alpha1 * I(covar_bin==1)+alpha2 * covar_conti + Aeps %*% matrix(betaA, ncol = 1) + Xeps %*% matrix(betaP, ncol = 1) + Meps %*% matrix(betaM, ncol = 1) + rnorm(1, mean = 0, sd = sqrt( 1 - (hsqA + hsqP + hsqM)))
    pheno <- as.integer(I (liability > thres_single ))
    
    if(pheno == 1 & (nrow(case_mat_ped_snps) - 1) < nobs*P ){
      case_mat_ped_snps <- rbind(case_mat_ped_snps, simulate_one_struct$mat_ped_snps)
      case_pat_ped_snps <- rbind(case_pat_ped_snps, simulate_one_struct$pat_ped_snps)
      case_mat_ped2_snps <- rbind(case_mat_ped2_snps, simulate_one_struct$mat_ped2_snps)
      tag <- tag + 1
      case_covar_bin <- c(case_covar_bin,covar_bin)
      case_covar_conti <- c(case_covar_conti,covar_conti)
    }
    if(pheno == 0 & (nrow(control_mat_ped_snps) - 1) < nobs*(1-P)){
      control_mat_ped_snps <- rbind(control_mat_ped_snps,simulate_one_struct$mat_ped_snps) 
      control_pat_ped_snps <- rbind(control_pat_ped_snps,simulate_one_struct$pat_ped_snps) 
      control_mat_ped2_snps <- rbind(control_mat_ped2_snps,simulate_one_struct$mat_ped2_snps) 
      tag <- tag + 1
      control_covar_bin <- c(control_covar_bin,covar_bin)
      control_covar_conti <- c(control_covar_conti,covar_conti)
      
    }

  }
  mat_ped_snps <- rbind(case_mat_ped_snps[-1,],control_mat_ped_snps[-1,])
  pat_ped_snps <- rbind(case_pat_ped_snps[-1,],control_pat_ped_snps[-1,])
  mat_ped2_snps <- rbind(case_mat_ped2_snps[-1,],control_mat_ped2_snps[-1,])
  
  simulate_struct <- list()
  
  simulate_struct$mat_ped_snps <- mat_ped_snps
  simulate_struct$pat_ped_snps <- pat_ped_snps
  simulate_struct$mat_ped2_snps <- mat_ped2_snps
  
  covar1 <- c(case_covar_bin,control_covar_bin)
  covar2 <- c(case_covar_conti,control_covar_conti)
  
  # generate the .ped and .map files in PLINK
  with(simulate_struct, make_mldose_and_ped(mat_ped_snps = mat_ped_snps, pat_ped_snps = pat_ped_snps, mat_ped2_snps = mat_ped2_snps,runid = times, nobs = nobs, nsnps = nsnps,dir=dir,simname=simname ))
  make_adp_grms(times,dir=dir,simname=simname)   
  
  pheno <- c(rep(1,nobs*P),rep(0,nobs*(1-P)))
  write.table(cbind(paste0("id", 1:nobs),  paste0("id", 1:nobs),  pheno), file = paste0(dir,simname,"_",times,".pheno"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(cbind(paste0("id", 1:nobs),  paste0("id", 1:nobs),  covar1), file = paste0(dir,simname,"_",times,".covar"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(cbind(paste0("id", 1:nobs),  paste0("id", 1:nobs),  covar2), file = paste0(dir,simname,"_",times,".qcovar"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  rm(list = ls())
  
  }
      


estima_CC <- function(dir1,times,nobs,K,P=0.5){
  library(quadprog)
  source("~/simulate_funcs_apm.R")
  
  dir=dir1
  simname <- paste0("bin_",K)
  
  # REML
  system(paste0("touch simgrmlist.", times, ".txt"))
  system(paste0("echo " ,dir, "simA.", times, ">> simgrmlist.", times, ".txt"))     
  system(paste0("echo ",dir,"simP.", times, ">> simgrmlist.", times, ".txt"))
  system(paste0("echo ",dir, "simM.", times, ".d >> simgrmlist.", times, ".txt"))
  system(paste0("mv simgrmlist.",times,".txt ",dir))
  
  system(paste0(cmd_gcta,
                paste0("./gcta64 --reml --reml-no-constrain  --mgrm-bin ",dir,"simgrmlist.", times,
                       ".txt --prevalence ", K," --pheno ",dir ,simname,"_",times,
                       ".pheno  --covar ",dir ,simname,"_",times,
                       ".covar --qcovar ",dir ,simname,"_",times,
                       ".qcovar --out ",dir,simname,"_APM_covar_", times)), 
         intern = TRUE)
  
  # PCGC & R-PCGC
  
  prefix=paste0(dir, "simA.", times)
  result=ReadGRMBin(prefix)
  A_c_n <- result$off

  prefix=paste0(dir, "simP.", times)
  result=ReadGRMBin(prefix)
  P_c_n <- result$off

  prefix=paste0(dir, "simM.", times,".d")
  result=ReadGRMBin(prefix)
  M_c_n <- result$off
  
  ### read the phenotype file and calculate the phenotypic correlations
  phen <- read.table(paste0(dir,simname,"_",times,".pheno"),header = F,sep=" ")[,3]
  covar1 <- read.table(paste0(dir,simname,"_",times,".covar"),header=F,sep=" ")[,3]
  covar2 <- read.table(paste0(dir,simname,"_",times,".qcovar"),header=F,sep=" ")[,3]
  C_c_multi <- transformPheno(pheno,cbind(factor(covar1),covar2),K)$multiplier
  w <- transformPheno(pheno,cbind(factor(covar1),covar2),K)$phen## the input outcome for the linear regression
  
  C_c_mat <- C_c_multi %*% t(C_c_multi)
  C_c_coef <- as.vector(C_c_mat[upper.tri(C_c_mat,diag = FALSE)])
  
  m_w <- w %*% t(w) ## the input outcome for the linear regression
  p_c <- as.vector(m_w[upper.tri(m_w,diag = FALSE)]) ## the phenotype correlation vector
  rm("m_w")
  
  A_c_fix <- A_c_n * C_c_coef
  P_c_fix <- P_c_n * C_c_coef
  M_c_fix <- M_c_n * C_c_coef
  
  rm("A_c_n","P_c_n","M_c_n","C_c_mat","C_c_coef")
  
  fit_pcgc_fix <- lm(p_c ~ A_c_fix + P_c_fix + M_c_fix)
  sink(paste0(dir,simname,"_PCGC_covar_", times,".txt"))
  print(summary(fit_pcgc_fix))
  sink()
  
  ## R-PCGC
  y_matrix = p_c
  X_matrix = cbind( A_c_fix , P_c_fix , M_c_fix)
  XTX = crossprod(X_matrix)
  XTY = crossprod(X_matrix, y_matrix)
  
  repcgc = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,3), bvec = rep(0, 3), meq=0, factorized=FALSE)$solution
  
  sum_XTX <- matrix(rowSums(apply(X_matrix,1,function(X) X%*%t(X))),nrow=3,byrow = F)
  inver_sum_XTX <- solve(sum_XTX)
  write.csv(c(repcgc,sqrt(diag(inver_sum_XTX ))),paste0(dir,simname,"_RE-PCGC_",times,"_covar.txt"))
  rm(list = ls())
}
 