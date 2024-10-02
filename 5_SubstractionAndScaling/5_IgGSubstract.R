################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 5. Normalization                                                             #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: substract igg remove negative values                                #
################################################################################

#-------------------Paths------------------------------------------------------#

args = commandArgs(trailingOnly=TRUE)
print(length(args))
root = args[1]
bwpath <- paste0(root, "/alignment/bigwig/")
file <- args[2]
dupremove <- args[3]
#----------------Igg substaction-----------------------------------------------#

library("GenomicRanges")
library("rtracklayer")

iggzeroes <- function(file, bwpath){
  cat(file)
  gr <- import(paste0(bwpath, file),
               as = "GRanges")
  
  gr$score[gr$score < 0] <- 0 
  
  export(gr, paste0(bwpath, file),
         format = "bw")
  return(0)
}

filenorm <- paste0(file, "_bowtie2.fragments.normalized.substracted.igg.bw")
iggzeroes(filenorm, bwpath)
filefrag <- paste0(file, "_bowtie2.fragments.substracted.igg.bw")
iggzeroes(filefrag, bwpath)
if (dupremove == "true"){
  filededup <- paste0(file, "_bowtie2.rmDup.fragments.substracted.igg.bw")
  iggzeroes(filededup, bwpath)
  filenormdedup <- paste0(file, "_bowtie2.rmDup.fragments.normalized.substracted.igg.bw")
  iggzeroes(filenormdedup, bwpath)
}

  

