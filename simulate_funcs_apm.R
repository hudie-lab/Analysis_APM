## janehoe@mail.ustc.edu.cn
## 06/25/2022

## Function to generate the GRMs for APM model on case-control studies 
## due to additive, parent-of-origin and maternal genetic effects.

cmd_plink <- "cd ~/plink;" ## the linux terminal command for open the file of PLINK
cmd_gcta <- "cd ~/gcta;" ## the linux terminal command for open the file of GCTA

## The functions follows the work of Laurin et al. (2018) (https://github.com/amatrhr/g-remladp)

make_mldose_and_ped <- function(nsnps, nobs, mat_ped_snps, pat_ped_snps,mat_ped2_snps, times,dir,simname) {
    
    sim_ids <- paste0("id", 1:nobs)
    sim_snps <- paste0("rs", 1:nsnps)
    
    ## MLDOSE
    ped_snps <- mat_ped_snps + pat_ped_snps
    dosage <- mat_ped_snps - pat_ped_snps
    
    ## MLINFO
    mlmeans <- colMeans(ped_snps)/2
    mlA1 <- rep("A", nsnps)
    mlA2 <- rep("G", nsnps)
    mlQ <- runif(n = nsnps, min = 0.5, max = 0.999) 
    mlRsq <- runif(n = nsnps, min = 0.5, max = 0.999)

    mlinfo <- data.frame(sim_snps, mlA1, mlA2, mlmeans, mlQ, mlRsq)

    write.table(data.frame(paste(sim_ids, sim_ids, sep = "->"),"ML_DOSE", dosage), file = gzfile(paste0(dir,simname, times, ".mldose.gz")), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(data.frame(sim_snps, mlA1, mlA2, sapply(1:nsnps, FUN = function(i) {min(mlmeans[i], 1- mlmeans[i])}), mlmeans, mlQ, mlRsq), file = gzfile(paste0(dir,simname, times, ".mlinfo.gz")), row.names = FALSE, col.names = c("SNP", "Al1", "Al2", "Freq1", "MAF", "Quality", "Rsq"), sep = "\t", quote = FALSE)

    ## MAP
    sim_map <- cbind(1, sim_snps, 0, 10000 + 1:nsnps)
    write.table(sim_map, file = paste0(dir,simname,times,".map"), row.names = FALSE, col.names = FALSE, quote = FALSE)

    ## PED
    obs_sexes <- sample(1:2,nobs,replace = TRUE)
    ped_leading_cols <- cbind(sim_ids, sim_ids, 0, 0, obs_sexes, -9)
    ped_snps_towrite <- apply(ped_snps, 2, gsub, pattern = 0, replacement = "A A")
    ped_snps_towrite <- apply(ped_snps_towrite, 2, gsub, pattern = 1, replacement = "A G")
    ped_snps_towrite <- apply(ped_snps_towrite, 2, gsub, pattern = 2, replacement = "G G")
    write.table(cbind(ped_leading_cols, ped_snps_towrite), file = paste0(dir,simname,times,".ped"), row.names = FALSE, col.names = FALSE, quote = FALSE)

    ## Convert to binary and remove big-n-sloppy plaintext
    system(paste0(cmd_plink,paste0("./plink --noweb --file ",dir,simname, times," --recode --make-bed --out ",dir,simname, times)))

    #### for ped of mom
    
    ped_snps <- mat_ped_snps + mat_ped2_snps

    ## MAP
    #sim_map <- cbind(1, sim_snps, 0, 10000 + 1:nsnps)
    write.table(sim_map, file = paste0(dir,simname,"M.",times,".map"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    ## PED
    ped_leading_cols <- cbind(sim_ids, sim_ids, 0, 0, obs_sexes, -9)
    ped_snps_towrite <- apply(ped_snps, 2, gsub, pattern = 0, replacement = "A A")
    ped_snps_towrite <- apply(ped_snps_towrite, 2, gsub, pattern = 1, replacement = "A G")
    ped_snps_towrite <- apply(ped_snps_towrite, 2, gsub, pattern = 2, replacement = "G G")
    write.table(cbind(ped_leading_cols, ped_snps_towrite), file = paste0(dir,simname,"M.",times,".ped"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    ## Convert to binary and remove big-n-sloppy plaintext
    system(paste0(cmd_plink,paste0("./plink --noweb --file ",dir,simname,"M.", times," --recode --make-bed --out ",dir,simname,"M.", times)))

}



make_adp_grms <- function(times,dir,simname){

  ## Additive GRM
  system(paste0(cmd_gcta,paste0("./gcta64 --bfile ",dir,simname, times, " --make-grm-bin --out ",dir,"simA.", times)))
  ## Parent-of-origin GRM
  system(paste0(cmd_gcta,paste0("./gcta64 --dosage-mach-gz ",dir ,simname, times, ".mldose.gz ",dir,simname, times, ".mlinfo.gz --make-grm-bin --out ",dir,"simP.",times)))
  ## Maternal GRM
    system(paste0(cmd_gcta,paste0("./gcta64 --bfile ",dir,simname,"M.", times, " --make-grm-d-bin --out ",dir,"simM.", times)))
  
}