################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Mapping Summary                                                           #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: R script to compute summary and statistics of the mappings          #
################################################################################

#-------------------------Paths------------------------------------------------#
library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]

#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(paste0(projPath, 
                                 "/experiment_summary_Latest.csv"),
                          header = T, sep = ",")

sampleList <- sampletable$SampleName

nameList <- unique(sapply(strsplit(sampleList,"_"), `[`, 1))

## Collect the alignment results from the bowtie2 alignment summary files
alignResult = c()
for(exp in sampleList){
  alignRes = read.table(paste0(projPath, 
                              "/alignment/sam/bowtie2_summary/", 
                              exp, "_bowtie2.txt"), 
                       header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  expInfo = strsplit(exp, "_")[[1]]
  print(expInfo)
  alignResult = data.frame(Experiment = expInfo[1], Replicate = expInfo[2], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% 
                                             as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% 
                                                as.character %>% as.numeric + 
                                                alignRes$V1[5] %>% as.character %>% 
                                                as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% 
                                                rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Experiment, levels = nameList)
alignResult <- alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, 
                                                           "%"))
write.table(alignResult, paste0(projPath,
                         "/alignment/summary_seq_depth_all_experiments.txt"),
                         row.names = FALSE)
