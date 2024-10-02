################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Dup Assess Summary                                                        #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: R script to compute summary and statistics of the duplicates        #
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
                           header=T)
alignSummary$Replicate <- as.character(alignSummary$Replicate)
sampleList <- sampletable$SampleName

nameList <- unique(sapply(strsplit(sampleList,"_"), `[`, 1))

## Collect the alignment results from the bowtie2 alignment summary files
dupResult = c()
for(exp in sampleList){
  dupRes = read.table(paste0(projPath, 
                             "/alignment/removeDuplicate/picard_summary/",
                             exp, "_picard.dupMark.txt"), 
                      header = TRUE, 
                      fill = TRUE)
  
  expInfo = strsplit(exp, "_")[[1]]
  dupResult = data.frame(Experiment = expInfo[1],
                         Replicate = expInfo[2], 
                         MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% 
                           as.character %>% as.numeric,
                         DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>%
                           as.character %>% as.numeric * 100,
                         EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>%
                           as.character %>% as.numeric) %>%
    mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% 
    rbind(dupResult, .)
}

dupResult$Experiment = factor(dupResult$Experiment, levels = nameList)
alignDupSummary = left_join(alignSummary, 
                            dupResult, 
                            by = c("Experiment", 
                                   "Replicate", 
                                   "MappedFragNum_hg38")) %>%
  mutate(DuplicationRate = paste0(DuplicationRate, "%"))


write.table(alignDupSummary, paste0(projPath,
                         "/alignment/summary_duplicates_all_experiments.txt"),
                         row.names = FALSE)


