################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 4. Reproducibility assessment                                                #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create table for reproducibility assessment plot                    #
################################################################################

library(dplyr, quietly=TRUE)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]
dupRemove <- args[3]
print(dupRemove)
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")

cat("Reproducibility of the replicates\n")

sampleList <- sort(sampletable$SampleName)
reprod = c()
fragCount = NULL
for(hist in sampleList){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
    
  }else{
    
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")

write.table(M, paste0(projPath,"/alignment/rep_reproducibility_all_experiments.txt"),
            row.names = T)

cat("Reproducibility of the replicates with duplicate removal\n")

if(dupRemove == "true"){
  sampleList <- sort(sampletable$SampleName)
  reprod = c()
  fragCount = NULL
  for(hist in sampleList){
    
    if(is.null(fragCount)){
      
      fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.rmDup.fragmentsCount.bin500.bed"), header = FALSE) 
      colnames(fragCount) = c("chrom", "bin", hist)
      
    }else{
      
      fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.rmDup.fragmentsCount.bin500.bed"), header = FALSE)
      colnames(fragCountTmp) = c("chrom", "bin", hist)
      fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
      
    }
  }
  
  M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")
  
  write.table(M, paste0(projPath,"/alignment/rep_reproducibility_all_experiments_rmDup.txt"),
              row.names = T)
  
}

