################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Mapping Summary Spike In                                                  #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: R script to compute summary and statistics of the mappings          #
################################################################################

#-------------------------Paths------------------------------------------------#
library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]


#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")

sampleList <- sampletable$SampleName

nameList <- unique(sapply(strsplit(sampleList,"_"), `[`, 1))


spikeAlign = c()
for(hist in sampleList){
  spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
  histInfo = strsplit(hist, "_", fixed = T)[[1]]
  spikeAlign = data.frame(Histone = histInfo[1], Replicate = paste(histInfo[-1], collapse = "_"), 
                          SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
                          MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
                          AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
}
spikeAlign$Histone = factor(spikeAlign$Histone, levels = nameList)
spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
write.table(spikeAlign, paste0(projPath,
                                "/alignment/summary_seq_depth_spikein.txt"),
            row.names = FALSE)


