
library(quadprog)

# input the directory of the GRM and the phenotype files and the prevalence K
dir <- "~/."

# function to read the GRM binary file
# ReadGRMBin is written in https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM

# R script to implement the methods of REML, PCGC and R-PCGC for case-control studies

## REML
system(paste0("touch grmlist.txt"))
system(paste0("echo " ,dir, "simA.>> grmlist.txt"))     
system(paste0("echo ",dir,"simP.>> grmlist.txt"))
system(paste0("echo ",dir, "simM.d >> grmlist.txt"))
system(paste0("mv grmlist.txt ",dir))

system(paste0(cmd_gcta,
              paste0("./gcta64 --reml --reml-no-constrain  --mgrm-bin ",dir,"grmlist.txt --pheno ",dir ,"_.pheno  --prevalence ", K , " --out ",dir,"_APM_")), 
       intern = TRUE)

## PCGC & R-PCGC

### read the GRMs and obetain the off-diagonal entries
prefix=paste0(dir, "simA.")
result=ReadGRMBin(prefix)
A_c_n <- result$off

prefix=paste0(dir, "simP.")
result=ReadGRMBin(prefix)
P_c_n <- result$off

prefix=paste0(dir, "simM.d")
result=ReadGRMBin(prefix)
M_c_n <- result$off

### read the phenotype file and calculate the phenotypic correlations
phen <- read.table(paste0(dir,".pheno"),header = F,sep=" ")[,3]
fit0  <- glm(phen ~ 1,family = "binomial")
w <- (phen-fit0$fitted.values)/sqrt(fit0$fitted.values*(1-fit0$fitted.values)) ## the input outcome for the linear regression
m_w <- w %*% t(w)
p_c <- as.vector(m_w[upper.tri(m_w,diag = FALSE)]) ## the phenotype correlation vector

### adjust the genetic correlations with the prevalence
P <- mean(phen)
adjust_coef = K^2*(1-K)^2/P/(1-P)/(dnorm(qnorm(1-K)))^2
A_c_pcgc <- A_c_n/adjust_coef
P_c_pcgc <- P_c_n/adjust_coef
M_c_pcgc <- M_c_n/adjust_coef


### PCGC 
fit <- lm(p_c ~ A_c_pcgc + P_c_pcgc + M_c_pcgc) ## Haseman-Elston regression
summary(fit) 
sink(paste0(dir,"_PCGC.txt"))
print(summary(fit))
sink()

### R-PCGC
y_matrix = p_c
X_matrix = cbind( A_c_pcgc , P_c_pcgc , M_c_pcgc)
XTX = crossprod(X_matrix)
XTY = crossprod(X_matrix, y_matrix)
repcgc = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,3), bvec = rep(0, 3), meq=0, factorized=FALSE)$solution

### the asymptic variance of R-PCGC estimators
sum_XTX <- matrix(rowSums(apply(X_matrix,1,function(X) X%*%t(X))),nrow=3,byrow = F)
inver_sum_XTX <- solve(sum_XTX)
write.table(c(repcgc,sqrt(diag(inver_sum_XTX ))),paste0(dir,"_R-PCGC_asymptotic.txt"),col.names = c("esti_A","esti_P","esti_M","sd_A","sd_P","sd_M"))
                                
                                
## to run analysis in parallel for 100 times
library(quadprog)

# input the directory of the GRM and the phenotype files and the prevalence K
N <- 3000
K <- 0.5
simname <- paste0("bin_", K)
dir <- "./"
#Cc_generate(dir,1, hsqA=0.4,hsqP = 0.05,hsqM = 0.1, K=K,nsnps = 5000, nobs=N)

sapply( 1:100, Cc_generate, dir=dir, hsqA=0.4,hsqP = 0.05,hsqM = 0.1, K=K,nsnps = 5000, nobs=N)

# estima_CC(dir1=dir, times=24, nobs=N,K=K)
sapply(1:100, estima_CC, dir1=dir, nobs=N, K=K)

## to read and restore the data
reml_result <- pcgc_result <- rpcgc_result <- matrix(nrow=108, ncol=6)
for(times in 1:100){
  varGCTA <- read.table(paste0(dir, simname,"_APM_", times, ".hsq"), header = F, sep = " ")[[1]]
  
  reml_result[times, 1:2] <- as.numeric(unlist(strsplit(varGCTA[(grepl("V.(G1.)/Vp_L", varGCTA))], "\t"))[2:3])
  reml_result[times, 3:4] <- as.numeric(unlist(strsplit(varGCTA[(grepl("V.(G2.)/Vp_L", varGCTA))], "\t"))[2:3])
  reml_result[times, 5:6] <- as.numeric(unlist(strsplit(varGCTA[(grepl("V.(G3.)/Vp_L", varGCTA))], "\t"))[2:3])
  
  varPCGC <- read.table(paste0(dir, simname, "_PCGC_", times, ".txt"), header = F, sep = ",")
  
  pcgc_result[times, 1:2] <- as.numeric(unlist(strsplit(varPCGC[9,], " "))[
    nchar(unlist(strsplit(varPCGC[9,], " "))) > 0 ][2:3])
  pcgc_result[times, 3:4] <- as.numeric(unlist(strsplit(varPCGC[10,], " "))[
    nchar(unlist(strsplit(varPCGC[10,], " "))) > 0 ][2:3])
  pcgc_result[times, 5:6] <- as.numeric(unlist(strsplit(varPCGC[11,], " "))[
    nchar(unlist(strsplit(varPCGC[11,], " "))) > 0 ][2:3])
  
  varRpcgc <- as.numeric(unlist(
    read.table(paste0(dir, simname, "_R-PCGC_", times, ".txt"), header = F, sep = ",")[,-1])[-1])
  
  rpcgc_result[times, 1:2] <- varRpcgc[1:2]
  rpcgc_result[times, 3:4] <- varRpcgc[3:4]
  rpcgc_result[times, 5:6] <- varRpcgc[5:6]
}

reml_result <- reml_result[-which(rowSums(is.na(reml_result))!=0),]
pcgc_result <- pcgc_result[-which(rowSums(is.na(pcgc_result))!=0),]
rpcgc_result <- rpcgc_result[-which(rowSums(is.na(rpcgc_result))!=0),]

write.csv(reml_result, "./reml_result.csv")
write.csv(pcgc_result, "./pcgc_result.csv")
write.csv(rpcgc_result, "./rpcgc_result.csv")
