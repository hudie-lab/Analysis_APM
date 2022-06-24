## janehoe@mail.ustc.edu.cn
## 06/25/2022

## R script to recode the genotypes of children due to parent-of-origin genetic effects, 
## and to generate the .mlinfo.gz and .mldosage.gz files.

# the directory to store the files
dir <- "~/"

# the .fam file contains the SNP data of mother-child duos
data.fam <- read.table("MotherChild.fam")
# the .freq file contains the information of SNP data of mother-child duos
data.freq <- read.table("snp.frq",stringsAsFactors = F,header=T)

snp_poo_imprint <- function(ped_split){
  
  # ped_split is the family' ped information of one pair of mother-chiod duo
  
  a <- ped_split[5] # the minor allele for the ith SNP
  p <- ped_split[6]
  mat_snp <- sum(grepl(a, c(ped_split[1:2])))
  child_snp <- sum(grepl(a, c(ped_split[3:4])))
  mat_ped <- integer() # the allele transparented from mother
  pat_ped <- integer() # the allele transparented from father
  
  if(mat_snp == 1){
    if(child_snp == 1){
      # the case when all the family members are heterozygous
      pat_ped <- rbinom(1,1,prob = p)
      mat_ped <- child_snp - pat_ped
    } 
    if(child_snp == 0){
      mat_ped <-0
    } 
    if(child_snp == 2){
      mat_ped <- 1
    }
    pat_ped <- child_snp - mat_ped
  }
  
  if(mat_snp == 0){
    mat_ped <- 0
    pat_ped <- child_snp - mat_ped
  } 
  if(mat_snp == 2 ){
    mat_ped <- 1
    pat_ped <- child_snp - mat_ped
  } 
  
  imprint <-  mat_ped - pat_ped # the mat_ped - pat_ped as the imprinting code
  return(imprint) 
}

# chromo_choice: the interested chromosomo to recode the genotypes
for(c in chromo_choice){
  write.table(data.freq[data.freq$CHR == c,"SNP"],paste0(dir,"snp_chromo_",c,".txt"),col.names = F,row.names = F,quote = F)
  freq_c <- data.freq[data.freq$CHR == c,]$A1
  freq_maf <- data.freq[data.freq$CHR == c,]$MAF
  system(paste0(plink_dir,"./plink --noweb --bfile MotherChild --extract ",dir,"snp_chromo_",c,".txt --make-bed --out ",dir,"single_chromo_samples_snp_",c))
  
  imprint_matrix <- data.frame()
  
  fam.list <- data.fam[,c("FID","IID")]
  for(f in fam.list$FID){
    
    family.list <- data.fam[data.fam$V1 == f, 1:2]
    write.table(family.list,paste0(dir,"single_fam_",f,".txt"),col.names = F,row.names = F,quote = F)
    
    system(paste0("./plink --noweb --bfile ",dir,"single_chromo_samples_snp_",c," --keep ",dir,"single_fam_",f,".txt --recode --out ",dir,"single_fam_samples_snp_",c))
    
    data.ped <-  read.table(paste0(dir,"single_fam_samples_snp_",c,".ped"),header = F,stringsAsFactors = F,colClasses = "character")
    fam_id <- data.ped[data.ped$V1 == f, 1:4] # the id of the family
    child_id <- fam_id[is.na(match(fam_id$V2,fam_id$V4)),] # the id of the child
    mat_id <- child_id[,c(1,4)] # the id of the mother
    ped_split <- rbind(matrix(unlist(data.ped[data.ped$V2 == mat_id[2][[1]],-(1:6)]),byrow = F,nrow=2), 
                       matrix(unlist(data.ped[data.ped$V2 == child_id[2][[1]],-(1:6)]),byrow = F,nrow=2),
                       freq_c,freq_maf) # select the family' ped information
    
    split_result <- unlist(apply(ped_split,2,snp_poo_imprint))
    
    imprint_matrix <- rbind(imprint_matrix,split_result)
    rm(data.ped)
    
  }
  write.table(cbind(fam.list,imprint_matrix),paste0(dir,"imprinting_code",c,".txt"),row.names = F,col.names = F,quote = F)
}

#----------------generate the dosage zip file------------------------------------------------# 
#--------------------------------------------------------------------------------------------#
dosage <- data.frame(fam.list)

for(c in unique(data.freq$CHR)){
  imprint_matrix <- read.table(paste0(dir,"imprinting_code",c,".txt"),header=F,stringsAsFactors = F)
  dosage <- cbind(dosage, imprint_matrix[,-(1:2)])
}
#--------------------------------------------------------------------------------------------#

#------------generate the dosage zip file----------------------------------------------------#
write.table(data.frame(paste(fam.list[,1], fam.list[,2], sep = "->"),"ML_DOSE", dosage[,-(1:2)]), file = gzfile(paste0(dir,"samples", ".mldose.gz")), row.names = FALSE, col.names = FALSE, quote = FALSE)
#--------------------------------------------------------------------------------------------#

#------------------generate the info zip file------------------------------------------------# 
## here the .mlinfo.gz is generated randomly
mlQ <- runif(n = nrow(data.freq), min = 0.5, max = 0.999) 
mlRsq <- runif(n = nrow(data.freq), min = 0.5, max = 0.999)
mlinfo <- data.frame(data.freq[,2:5], 1- data.freq$MAF, mlQ, mlRsq)
colnames(mlinfo) <- c("SNP", "Al1", "Al2", "Freq1", "MAF", "Quality", "Rsq")
write.table(mlinfo, file = gzfile(paste0(dir,"samples", ".mlinfo.gz")), row.names = FALSE, col.names = T, sep = "\t", quote = FALSE)
