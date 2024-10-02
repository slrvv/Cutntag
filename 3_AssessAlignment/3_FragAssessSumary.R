################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Frag Assess Summary                                                       #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: R script to compute summary and statistics of the fragment sizes    #
################################################################################

#-------------------------Paths------------------------------------------------#
library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]

#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(paste0(summaryPath),
                          header = T, sep = ",")
alignSummary <- read.table(paste0(projPath,
                                  "/alignment/summary_seq_depth_all_experiments.txt"),
                           header=T, sep=",")
sampleList <- sampletable$SampleName

nameList <- unique(sapply(strsplit(sampleList,"_"), `[`, 1))

fragLen = c()

for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[1], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}

fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = nameList)
write.table(fragLen, paste0(projPath,
                            "/alignment/summary_fraglen_all_experiments.txt"),
            row.names = FALSE)



