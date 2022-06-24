## janehoe@mail.ustc.edu.cn
## 06/25/2022

## function to obtain the REML estimators on case-control studies of the heritabilities 
## due to additive, parent-of-origin and maternal genetic effects from the software GCTA.

REML <- function(dir,K){
  # dir: the directory of files
  # K: the prevalence

  system(paste0("touch grmlist.txt"))
  system(paste0("echo " ,dir, "simA.>> grmlist.txt"))     
  system(paste0("echo ",dir,"simP.>> grmlist.txt"))
  system(paste0("echo ",dir, "simM.d >> grmlist.txt"))
  system(paste0("mv grmlist.txt ",dir))
  
  system(paste0(cmd_gcta,
                paste0("./gcta64 --reml --reml-no-constrain  --mgrm-bin ",
                       dir,"grmlist.txt --prevalence ", K," --pheno ",
                       dir ,".pheno  --covar ",dir ,
                       ".covar --qcovar ",dir ,
                       ".qcovar --out ",dir,"_APM_covar")), 
         intern = TRUE)

}