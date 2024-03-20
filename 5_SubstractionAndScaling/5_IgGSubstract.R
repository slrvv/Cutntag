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
Bedgraphpath <- paste0(root, "/alignment/bigwig/")
FILE1 <- args[2]

#----------------Igg substaction-----------------------------------------------#

library("GenomicRanges")
library("rtracklayer")

#gr <- import(paste0(Bedgraphpath, FILE1), as = "GRanges")
print(FILE1)
print(paste0(Bedgraphpath, FILE1))
gr <- import(paste0(Bedgraphpath, FILE1),
             as = "GRanges")

gr$score[gr$score < 0] <- 0 

export(gr, paste0(Bedgraphpath, FILE1),
             format = "bw")
