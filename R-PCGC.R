## janehoe@mail.ustc.edu.cn
## 06/25/2022

## function to obtain the R-PCGC estimators on case-control studies of the heritabilities 
## due to additive, parent-of-origin and maternal genetic effects.

R-PCGC <- function(dir,K){
  # dir: the directory of files
  # K: the prevalence

  # ReadGRMBin() is written in https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
  ### read the GRMs and obetain the off-diagonal entries
  prefix=paste0(dir, "simA.", times)
  result=ReadGRMBin(prefix)
  A_c_n <- result$off
  
  prefix=paste0(dir, "simP.", times)
  result=ReadGRMBin(prefix)
  P_c_n <- result$off
  
  prefix=paste0(dir, "simM.", times,".d")
  result=ReadGRMBin(prefix)
  M_c_n <- result$off
  
  # transformPheno() is written by David Golan in https://github.com/gauravbhatia1/PCGCRegression
  ### read the phenotype file and calculate the phenotypic correlations
  phen <- read.table(paste0(dir,simname,"_",times,".pheno"),header = F,sep=" ")[,3]
  covar <- read.table(paste0(dir,simname,"_",times,".covar"),header=F,sep=" ")[,3]
  C_c_multi <- transformPheno(pheno,covar,K)$multiplier
  w <- transformPheno(pheno,covar,K)$phen## the input outcome for the linear regression
  
  C_c_mat <- C_c_multi %*% t(C_c_multi)
  C_c_coef <- as.vector(C_c_mat[upper.tri(C_c_mat,diag = FALSE)])
  
  m_w <- w %*% t(w) ## the input outcome for the linear regression
  p_c <- as.vector(m_w[upper.tri(m_w,diag = FALSE)]) ## the phenotype correlation vector
  rm("m_w")
  
  A_c_fix <- A_c_n * C_c_coef
  P_c_fix <- P_c_n * C_c_coef
  M_c_fix <- M_c_n * C_c_coef
  
  rm("A_c_n","P_c_n","M_c_n","C_c_mat","C_c_coef")
 
  ## R-PCGC
  y_matrix = p_c
  X_matrix = cbind( A_c_fix , P_c_fix , M_c_fix)
  XTX = crossprod(X_matrix)
  XTY = crossprod(X_matrix, y_matrix)
  repcgc = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,3), bvec = rep(0, 3), meq=0, factorized=FALSE)$solution
  
  ### to estimate asymptic variance of R-PCGC estimators
  sum_XTX <- matrix(rowSums(apply(X_matrix,1,function(X) X%*%t(X))),nrow=3,byrow = F)
  inver_sum_XTX <- solve(sum_XTX)
  
  return(list(est_h2=repcgc, esti_asyp_sd=sqrt(diag(inver_sum_XTX )) ))

}